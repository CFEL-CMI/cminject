/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause
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

/* forcedPoiseuille2d.cpp:
 * In this implementation of a Poiseuille flow, the boundaries are
 * periodic between inlet and outlet. The flow is driven by a body
 * force.
 * It illustrates the use of a body force,
 * periodic boundaries and computation of error norms.
 */


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR ForcedD2Q9Descriptor


// Parameters for the simulation setup
const T lx  = 2.;       // length of the channel
const T ly  = 1.;       // height of the channel
const int N = 50;       // resolution of the model
const int M = N;        // time discretization refinement
const T Re = 10.;       // Reynolds number
const T maxPhysT = 20.; // max. simulation time in s, SI unit


// Stores geometry information in form of material numbers
void prepareGeometry( LBconverter<T> const& converter,
                      SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );

  Vector<T,2> extend( lx, ly - 1.8*converter.getLatticeL() );
  Vector<T,2> origin( T(), 0.9*converter.getLatticeL() );
  IndicatorCuboid2D<T> cuboid2( extend, origin );

  superGeometry.rename( 2,1,cuboid2 );

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
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  T   omega = converter.getOmega();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 2, &bulkDynamics );

  // Setting of the boundary conditions
  sBoundaryCondition.addVelocityBoundary( superGeometry, 2, omega );

  // Initial conditions
  T Ly = converter.numCells( ly );
  std::vector<T> poiseuilleForce( 2,T() );
  poiseuilleForce[0] = 8.*converter.getLatticeNu()*converter.getLatticeU() / ( Ly*Ly );
  AnalyticalConst2D<T,T> force( poiseuilleForce );

  // Initialize force
  sLattice.defineExternalField( superGeometry, 1,
                                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                DESCRIPTOR<T>::ExternalField::sizeOfForce, force );
  sLattice.defineExternalField( superGeometry, 2,
                                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                DESCRIPTOR<T>::ExternalField::sizeOfForce, force );

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
  const T radius = ly/2.;
  std::vector<T> axisPoint( 2,T() );
  axisPoint[0] = lx/2.;
  axisPoint[1] = ly/2.;
  std::vector<T> axisDirection( 2,T() );
  axisDirection[0] = 1;
  axisDirection[1] = 0;
  Poiseuille2D<T> uSol( axisPoint, axisDirection, maxVelocity, radius );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> uSolLattice( uSol,sLattice,superGeometry );

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

  // strainRate error
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

  SuperVTMwriter2D<T> vtmWriter( "forcedPoiseuille2d" );
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
      const T radius = ly/2.;
      std::vector<T> axisPoint( 2,T() );
      axisPoint[0] = lx/2.;
      axisPoint[1] = ly/2.;
      std::vector<T> axisDirection( 2,T() );
      axisDirection[0] = 1;
      axisDirection[1] = 0;
      Poiseuille2D<T> uSol( axisPoint, axisDirection, maxVelocity, radius );
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

  LBconverter<T> converter(
    ( int ) 2,                             // dim
    1./N,                                  // latticeL_
    1./M,                                  // latticeU_
    ( T )   1./Re,                         // charNu_
    ( T )   1.                             // charL_ = 1,
  );
  converter.print();
  // Writes the converter log file
  writeLogFile( converter, "forcedPoiseuille2d" );

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

  // Periodic boundaries in x-direction
  cuboidGeometry.setPeriodicity( true, false );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  ForcedBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getOmega(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition2D<T, DESCRIPTOR> ( sBoundaryCondition );
  //createLocalBoundaryCondition2D<T, DESCRIPTOR> (sBoundaryCondition);

  prepareLattice( converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry );

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

