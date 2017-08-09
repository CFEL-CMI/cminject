/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* bgkPoiseuille2d.cpp:
 * This example examines a 2D Poseuille flow with a velocity
 * or pressure boundary at the inlet/outlet.
 * Computation of error norms via functors is shown as well.
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef enum {velocity, pressure} BoundaryType;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor


// Parameters for the simulation setup
const T lx  = 2.;      // length of the channel
const T ly  = 1.;      // height of the channel
const int N = 100;      // resolution of the model
const int M = N;       // time discretization refinement
const T Re = 10.;      // Reynolds number
const T maxPhysT = 5.; // max. simulation time in s, SI unit


// Stores geometry information in form of material numbers
void prepareGeometry( LBconverter<T> const& converter,
                      SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,1,1 );

  Vector<T,2> extend;
  Vector<T,2> origin;

  // Set material number for inflow
  extend[1] = ly;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = lx;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( LBconverter<T> const& converter,
                     SuperLattice2D<T, DESCRIPTOR>& sLattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                     BoundaryType inflowBoundary, BoundaryType outflowBoundary,
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  T   omega = converter.getOmega();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 2, &bulkDynamics );

  // Material=3 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );

  // Material=4 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );

  // Setting of the boundary conditions
  sBoundaryCondition.addVelocityBoundary( superGeometry, 2, omega );

  if ( inflowBoundary==velocity ) {
    sBoundaryCondition.addVelocityBoundary( superGeometry, 3, omega );
  } else {
    sBoundaryCondition.addPressureBoundary( superGeometry, 3, omega );
  }

  if ( outflowBoundary==velocity ) {
    sBoundaryCondition.addVelocityBoundary( superGeometry, 4, omega );
  } else {
    sBoundaryCondition.addPressureBoundary( superGeometry, 4, omega );
  }

  // Initial conditions
  T Lx = converter.numCells( lx );
  T Ly = converter.numCells( ly );

  T p0 =8.*converter.getLatticeNu()*converter.getLatticeU()*Lx/( Ly*Ly );
  AnalyticalLinear2D<T,T> rho( -p0/lx*DESCRIPTOR<T>::invCs2 , 0 , p0*DESCRIPTOR<T>::invCs2+1 );

  T maxVelocity = converter.getLatticeU();
  T distance2Wall = converter.getLatticeL();
  Poiseuille2D<T> u( superGeometry, 3, maxVelocity, distance2Wall );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( superGeometry, 1, rho, u );
  sLattice.iniEquilibrium( superGeometry, 1, rho, u );
  sLattice.defineRhoU( superGeometry, 2, rho, u );
  sLattice.iniEquilibrium( superGeometry, 2, rho, u );
  sLattice.defineRhoU( superGeometry, 3, rho, u );
  sLattice.iniEquilibrium( superGeometry, 3, rho, u );
  sLattice.defineRhoU( superGeometry, 4, rho, u );
  sLattice.iniEquilibrium( superGeometry, 4, rho, u );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
void error( SuperGeometry2D<T>& superGeometry,
            SuperLattice2D<T, DESCRIPTOR>& sLattice,
            LBconverter<T> const& converter,
            Dynamics<T, DESCRIPTOR>& bulkDynamics ) {

  OstreamManager clout( std::cout,"error" );

  int tmp[]= {int()};
  T result[2]= {T(),T()}, result_tmp[2]= {T(),T()};
  T result1;

  const T maxVelocity = converter.physVelocity( converter.getLatticeU() );
  T distance2Wall = converter.getLatticeL();
  Poiseuille2D<T> uSol( superGeometry, 3, maxVelocity, distance2Wall );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> uSolLattice( uSol,sLattice,superGeometry );

  T Lx = converter.numCells( lx );
  T Ly = converter.numCells( ly );
  T p0 =8.*converter.getLatticeNu()*converter.getLatticeU()*Lx/( Ly*Ly );
  AnalyticalLinear2D<T,T> pressureSol( -converter.physPressureFromRho( p0*DESCRIPTOR<T>::invCs2+1 )/lx , 0 , converter.physPressureFromRho( p0*DESCRIPTOR<T>::invCs2+1 ) );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice,converter );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> pressureSolLattice( pressureSol,sLattice,superGeometry );

  PoiseuilleStrainRate2D<T,T> sSol( converter, ly );
  SuperLatticePhysStrainRate2D<T,DESCRIPTOR> s( sLattice,converter );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> sSolLattice( sSol,sLattice,superGeometry );


  // velocity error
  SuperL1Norm2D<T> uL1Norm( uSolLattice-u,superGeometry,1 );
  SuperL1Norm2D<T> uSolL1Norm( uSolLattice,superGeometry,1 );
  uL1Norm( result,tmp );
  uSolL1Norm( result_tmp,tmp );
  result1=result[0]/result_tmp[0];
  clout << "velocity-L1-error(abs)=" << result[0] << "; velocity-L1-error(rel)=" << result1 << std::endl;

  SuperL2Norm2D<T> uL2Norm( uSolLattice-u,superGeometry,1 );
  SuperL2Norm2D<T> uSolL2Norm( uSolLattice,superGeometry,1 );
  uL2Norm( result,tmp );
  uSolL2Norm( result_tmp,tmp );
  result1=result[0]/result_tmp[0];
  clout << "velocity-L2-error(abs)=" << result[0] << "; velocity-L2-error(rel)=" << result1 << std::endl;

  SuperLinfNorm2D<T> uLinfNorm( uSolLattice-u,superGeometry,1 );
  uLinfNorm( &result1,tmp );
  clout << "velocity-Linf-error(abs)=" << result1 << std::endl;

  // density error
  SuperL1Norm2D<T> pressureL1Norm( pressureSolLattice-pressure,superGeometry,1 );
  SuperL1Norm2D<T> pressureSolL1Norm( pressureSolLattice,superGeometry,1 );
  pressureL1Norm( result,tmp );
  pressureSolL1Norm( result_tmp,tmp );
  result1=result[0]/result_tmp[0];
  clout << "pressure-L1-error(abs)=" << result[0] << "; pressure-L1-error(rel)=" << result1 << std::endl;

  SuperL2Norm2D<T> pressureL2Norm( pressureSolLattice-pressure,superGeometry,1 );
  SuperL2Norm2D<T> pressureSolL2Norm( pressureSolLattice,superGeometry,1 );
  pressureL2Norm( result,tmp );
  pressureSolL2Norm( result_tmp,tmp );
  result1=result[0]/result_tmp[0];
  clout << "pressure-L2-error(abs)=" << result[0] << "; pressure-L2-error(rel)=" << result1 << std::endl;

  SuperLinfNorm2D<T> pressureLinfNorm( pressureSolLattice-pressure,superGeometry,1 );
  pressureLinfNorm( &result1,tmp );
  clout << "pressure-Linf-error(abs)=" << result1 << std::endl;

  // strain rate error
  SuperL1Norm2D<T> sL1Norm( sSolLattice-s,superGeometry,1 );
  SuperL1Norm2D<T> sSolL1Norm( sSolLattice,superGeometry,1 );
  sL1Norm( result,tmp );
  sSolL1Norm( result_tmp,tmp );
  result1=result[0]/result_tmp[0];
  clout << "strainRate-L1-error(abs)=" << result[0] << "; strainRate-L1-error(rel)=" << result1 << std::endl;

  SuperL2Norm2D<T> sL2Norm( sSolLattice-s,superGeometry,1 );
  SuperL2Norm2D<T> sSolL2Norm( sSolLattice,superGeometry,1 );
  sL2Norm( result,tmp );
  sSolL2Norm( result_tmp,tmp );
  result1=result[0]/result_tmp[0];
  clout << "strainRate-L2-error(abs)=" << result[0] << "; strainRate-L2-error(rel)=" << result1 << std::endl;

  SuperLinfNorm2D<T> sLinfNorm( sSolLattice-s,superGeometry,1 );
  sLinfNorm( &result1,tmp );
  clout << "strainRate-Linf-error(abs)=" << result1 << std::endl;
}

// Output to console and files
void getResults( SuperLattice2D<T,DESCRIPTOR>& sLattice, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 LBconverter<T>& converter, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer ) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "bgkPoiseuille2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int vtmIter  = converter.numTimeSteps( maxPhysT/20. );
  const int statIter = converter.numTimeSteps( maxPhysT/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    superGeometry.rename( 0,2 );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 ) {
    vtmWriter.write( iT );

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
    BlockLatticeReduction2D<T, DESCRIPTOR> planeReduction( normVel );
    BlockGifWriter<T> gifWriter;
    gifWriter.write( planeReduction, iT, "vel" ); // scaled

    ofstream *ofile = 0;
    if ( singleton::mpi().isMainProcessor() ) {
      ofile = new ofstream( ( singleton::directories().getLogOutDir()+"centerVel.dat" ).c_str() );
    }
    T Ly = converter.numCells( ly );
    for ( int iY=0; iY<=Ly; ++iY ) {
      T dx = converter.getDeltaX();
      const T maxVelocity = converter.physVelocity( converter.getLatticeU() );
      T point[2]= {T(),T()};
      point[0] = lx/2.;
      point[1] = ( T )iY/Ly;
      T distance2Wall = converter.getLatticeL();
      Poiseuille2D<T> uSol( superGeometry, 3, maxVelocity, distance2Wall );
      T analytical[2] = {T(),T()};
      uSol( analytical,point );
      SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
      AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> intpolateVelocity( velocity, true );
      T numerical[2] = {T(),T()};
      intpolateVelocity( numerical,point );
      if ( singleton::mpi().isMainProcessor() ) {
        *ofile << iY*dx << " " << analytical[0]
               << " " << numerical[0] << "\n";
      }
    }
    delete ofile;
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.physTime( iT ) );

    // Error norms
    error( superGeometry, sLattice, converter, bulkDynamics );
  }
}

int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  const BoundaryType inflowBoundary = velocity;
  const BoundaryType outflowBoundary = pressure;

  LBconverter<T> converter(
    ( int ) 2,                             // dim
    1./N,                                  // latticeL_
    1./M,                                  // latticeU_
    ( T )   1./Re,                         // charNu_
    ( T )   1.                             // charL_ = 1,
  );
  converter.print();
  writeLogFile( converter, "bgkPoiseuille2d" );

  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( lx, ly );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getLatticeL(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getOmega(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition2D<T,DESCRIPTOR>( sBoundaryCondition );
  // createLocalBoundaryCondition2D<T,DESCRIPTOR>(sBoundaryCondition);

  prepareLattice( converter, sLattice, bulkDynamics, sBoundaryCondition, inflowBoundary, outflowBoundary, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( int iT = 0; iT < converter.numTimeSteps( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, bulkDynamics, converter, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}

