/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2016 Vojtech Cvrcek, Mathias J. Krause
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

/* powerLaw2d.cpp:
 * This example examines a steady flow of a non-newtonian fluid in a channel.
 * At the inlet, a profile for non-newtonian fluid is imposed on the velocity,
 * where as the outlet implements an outflow condition grad_x p = 0.
 * The power law model is for n=1 and m=charNu in fact the classical poiseuille flow.
 * One can validate the error with using functors in void error.
 *
 *
 */

#include "olb2D.h"
#include "olb2D.hh"   // include full template code

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR DynOmegaD2Q9Descriptor

// Parameters for the simulation setup
int N = 40;            // resolution of the model
const T Re = 10.;      // Reynolds number
T lx = 2.;             // channel lenght
T ly = 1.;             // channel width
T maxPhysT = 20.;      // max. phys. time in s
// set the changes for n and m in powerLawBGKdynamics.h
T n = .2;              // parameter in power law model (n=1 Newtonian fluid)
T m = 1./Re;

void prepareGeometry( LBconverter<T> const& converter,
                      SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend( lx, ly );
  Vector<T,2> origin;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,1,1 );

  // Set material number for inflow
  extend[0] = 2.*converter.getLatticeL();
  origin[0] = -converter.getLatticeL();
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );
  // Set material number for outflow
  origin[0] = lx-converter.getLatticeL();
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice2D<T,DESCRIPTOR>& sLattice,
                     LBconverter<T> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getOmega();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  // Material=3 -->bulk dynamics (inflow)
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );

  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );

  // Setting of the boundary conditions
  sBoundaryCondition.addVelocityBoundary( superGeometry, 3, omega );
  sBoundaryCondition.addPressureBoundary( superGeometry, 4, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
}


void setBoundaryValues( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                        LBconverter<T> const& converter,
                        int iT, SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // Set initial and steady boundary conditions
  if ( iT==0 ) {

    // Define the analytical solutions for pressure and velocity
    T maxVelocity = converter.getLatticeU();
    T distance2Wall = converter.getLatticeL()/2.;

    T p0 = m*pow( converter.getCharU(),n )*pow( ( n + 1. )/n,n )*pow( 2./( ly-distance2Wall*2 ),n + 1. );

    AnalyticalLinear2D<T,T> rho( converter.rhoFromPhysicalPressure( -p0 ) - 1., 0, converter.rhoFromPhysicalPressure( p0*lx ) );

    PowerLaw2D<T> u( superGeometry, 3, maxVelocity, distance2Wall, ( n + 1. )/n );

    // Set the analytical solutions for pressure and velocity TODO???
    AnalyticalConst2D<T,T> omega0( converter.getOmega() );
    sLattice.defineExternalField ( superGeometry, 1,  DESCRIPTOR<T>::ExternalField::omegaBeginsAt, 1, omega0 );
    sLattice.defineExternalField ( superGeometry, 3, DESCRIPTOR<T>::ExternalField::omegaBeginsAt, 1, omega0 );
    sLattice.defineExternalField ( superGeometry, 4, DESCRIPTOR<T>::ExternalField::omegaBeginsAt, 1, omega0 );

    // Set the analytical solutions for pressure and velocity
    // Initialize all values of distribution functions to their local equilibrium

    sLattice.defineRhoU( superGeometry, 1, rho, u );
    sLattice.iniEquilibrium( superGeometry, 1, rho, u );

    sLattice.iniEquilibrium( superGeometry, 3, rho, u );
    sLattice.defineRhoU( superGeometry, 3, rho, u );

    sLattice.iniEquilibrium( superGeometry, 4, rho, u );
    sLattice.defineRhoU( superGeometry, 4, rho, u );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

// Compute error norms
void error( SuperGeometry2D<T>& superGeometry,
            SuperLattice2D<T, DESCRIPTOR>& sLattice,
            LBconverter<T> const& converter,
            Dynamics<T, DESCRIPTOR>& bulkDynamics ) {
  OstreamManager clout( std::cout,"error" );

  int tmp[1];
  T result[2];
  T result1[2];

  T distance2Wall = converter.getLatticeL()/2.;

  PowerLaw2D<T> uSol( superGeometry,3,converter.getCharU(),distance2Wall,( n + 1. )/n );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> uSolLattice( uSol,sLattice,superGeometry );

  T p0 = m*pow( converter.getCharU(),n )*pow( ( n + 1. )/n,n )*pow( 2./( ly-distance2Wall*2 ),n + 1. );
  AnalyticalLinear2D<T,T> pressureSol( -p0, 0, p0*lx );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> p( sLattice,converter );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> pSolLattice( pressureSol,sLattice,superGeometry );

  // velocity error
  SuperL1Norm2D<T> uL1Norm( uSolLattice-u,superGeometry,1 );
  SuperL1Norm2D<T> uSolL1Norm( uSolLattice,superGeometry,1 );
  uL1Norm( result,tmp );
  uSolL1Norm( result1,tmp );
  clout << "velocity-L1-error(abs)=" << result[0] << "; velocity-L1-error(rel)=" << result[0]/result1[0] << std::endl;

  SuperL2Norm2D<T> uL2Norm( uSolLattice-u,superGeometry,1 );
  SuperL2Norm2D<T> uSolL2Norm( uSolLattice,superGeometry,1 );
  uL2Norm( result,tmp );
  uSolL2Norm( result1,tmp );
  clout << "velocity-L2-error(abs)=" << result[0] << "; velocity-L2-error(rel)=" << result[0]/result1[0] << std::endl;

  SuperLinfNorm2D<T> uLinfNorm( uSolLattice-u,superGeometry,1 );
  SuperLinfNorm2D<T> uSolLinfNorm( uSolLattice,superGeometry,1 );
  uLinfNorm( result,tmp );
  uSolLinfNorm( result1,tmp );
  clout << "velocity-Linf-error(abs)=" << result[0] << "; velocity-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

  // pressure error
  SuperL1Norm2D<T> pL1Norm( pSolLattice-p,superGeometry,1 );
  SuperL1Norm2D<T> pSolL1Norm( pSolLattice,superGeometry,1 );
  pL1Norm( result,tmp );
  pSolL1Norm( result1,tmp );
  clout << "pressure-L1-error(abs)=" << result[0] << "; pressure-L1-error(rel)=" << result[0]/result1[0] << std::endl;

  SuperL2Norm2D<T> pL2Norm( pSolLattice-p,superGeometry,1 );
  SuperL2Norm2D<T> pSolL2Norm( pSolLattice,superGeometry,1 );
  pL2Norm( result,tmp );
  pSolL2Norm( result1,tmp );
  clout << "pressure-L2-error(abs)=" << result[0] << "; pressure-L2-error(rel)=" << result[0]/result1[0] << std::endl;

  SuperLinfNorm2D<T> pLinfNorm( pSolLattice-p,superGeometry,1 );
  SuperLinfNorm2D<T> pSolLinfNorm( pSolLattice,superGeometry,1 );
  pLinfNorm( result,tmp );
  pSolLinfNorm( result1,tmp );
  clout << "pressure-Linf-error(abs)=" << result[0] << "; pressure-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

}

// Output to console and files
void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                 Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 LBconverter<T> const& converter, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<double>& timer ) {
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "powerLaw2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  if ( iT==0 ) {
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice,superGeometry );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  if ( iT%converter.numTimeSteps( maxPhysT/20. )==0 ) {
    vtmWriter.write( iT );

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
    BlockLatticeReduction2D<T, DESCRIPTOR> planeReduction( normVel );
    BlockGifWriter<T> gifWriter;
    gifWriter.write( planeReduction, iT, "vel" ); // scaled
  }

  // Writes output on the console
  if ( iT%converter.numTimeSteps( maxPhysT/20. )==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT,iT*converter.getDeltaT() );
    error( superGeometry, sLattice, converter, bulkDynamics );
  }
  return;
}


int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  if ( argc > 1 ) {
    N = atoi( argv[1] );
  }
  if ( argc > 2 ) {
    n = atof( argv[2] );
  }

  singleton::directories().setOutputDir( "./tmp/" );

  LBconverter<T> converter(
    ( int ) 2,                      // dim
    ( T )   1./N,                   // latticeL_
    ( T )   1./N,                   // latticeU 0.05
    ( T )   m,                    // charNu_
    ( T )   1.,             // charL_ = 1
    ( T )   1.        // charU
  );
  converter.print();
  writeLogFile( converter, "powerLaw2d" );

  // === 2rd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights

  Vector<T,2> extend( lx, ly );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getLatticeL(), singleton::mpi().getSize() );
#else
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getLatticeL(), 7 );
#endif

  cuboidGeometry.printExtended();

  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );
  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  PowerLawBGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>(), m, n, converter.physTime() );

  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createLocalBoundaryCondition2D<T, DESCRIPTOR, PowerLawBGKdynamics<T,DESCRIPTOR> > ( sBoundaryCondition );

  prepareLattice( sLattice, converter, bulkDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  Timer<double> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( int iT=0; iT<converter.numTimeSteps( maxPhysT ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, bulkDynamics, converter, iT, superGeometry, timer );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
  }
  timer.stop();
  timer.printSummary();
}

