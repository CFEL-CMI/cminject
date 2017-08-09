/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

/* cavity3d.cpp:
 * This example illustrates a flow in a cuboid, lid-driven cavity.
 * This version is for parallel use. A version for sequential use
 * is also available.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif

#include <cmath>
#include <iostream>
#include <fstream>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor

const int N = 1; // resolution of the model
const int M = 1; // time discretization refinement
const T maxT = (T) 100.; // max. simulation time in s, SI unit

const T interval = 1.0; // Time intervall in seconds for convergence check
const T epsilon = 1e-3; // Residuum for convergence check

void prepareGeometry( LBconverter<T> const& converter, IndicatorF3D<T>& indicator, SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,1,1,1 );

  T eps = converter.getLatticeL();
  Vector<T,3> origin( -eps, converter.getCharL() - eps, -eps );
  Vector<T,3> extend( converter.getCharL() + 2*eps, 2*eps, converter.getCharL() + 2*eps );
  IndicatorCuboid3D<T> lid( extend,origin );

  superGeometry.rename( 2,3,1,lid );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( LBconverter<T> const& converter,
                     SuperLattice3D<T,DESCRIPTOR>& lattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getOmega();

  // Material=0 -->do nothing
  lattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  lattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2,3 -->bulk dynamics, velocity boundary
  lattice.defineDynamics( superGeometry, 2, &bulkDynamics );
  lattice.defineDynamics( superGeometry, 3, &bulkDynamics );
  bc.addVelocityBoundary( superGeometry, 2, omega );
  bc.addVelocityBoundary( superGeometry, 3, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( LBconverter<T> const&converter,
                        SuperLattice3D<T,DESCRIPTOR>& lattice, SuperGeometry3D<T>& superGeometry, int iT ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  if ( iT==0 ) {

    AnalyticalConst3D<T,T> rhoF( T( 1 ) );
    AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) );

    lattice.iniEquilibrium( superGeometry, 1, rhoF, uF );
    lattice.iniEquilibrium( superGeometry, 2, rhoF, uF );
    lattice.iniEquilibrium( superGeometry, 3, rhoF, uF );

    lattice.defineRhoU( superGeometry, 1, rhoF, uF );
    lattice.defineRhoU( superGeometry, 2, rhoF, uF );
    lattice.defineRhoU( superGeometry, 3, rhoF, uF );

    AnalyticalConst3D<T,T> uTop( converter.getLatticeU(), T( 0 ), T( 0 ) );
    lattice.defineU( superGeometry,3,uTop );

    // Make the lattice ready for simulation
    lattice.initialize();
  }
}

void getResults( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                 LBconverter<T> const& converter, SuperGeometry3D<T>& superGeometry,
                 int iT, Timer<T>& timer, bool converged ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "cavity3d" );

  const T logT     = ( T )1.;
  const T vtkSave  = ( T )1.;

  if ( iT==0 ) {
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( (iT%converter.numTimeSteps( logT )==0 && iT>0) || converged ) {
    timer.update( iT );
    timer.printStep( 2 );
    sLattice.getStatistics().print( iT,converter.physTime( iT ) );
  }

  // Writes the VTK and GIF files
  if ( (iT%converter.numTimeSteps( vtkSave )==0 && iT>0) || converged ) {
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    vtmWriter.write( iT );

    // define vector which span the gif-plane
    Vector<T,3> u( 1,0,0 );
    Vector<T,3> v( 0,1,0 );
    T tmp = T( converter.getCharL() / 2. );
    T origin[3] = {tmp,tmp,tmp};

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockLatticeReduction3D<T, DESCRIPTOR> planeReduction( normVel, u, v, 600, origin );
    BlockGifWriter<T> gifWriter;

    gifWriter.write( planeReduction, iT, "velocity" );
  }

  // Output for x-velocity along y-position at the last time step
  if ( (iT == converter.numTimeSteps( maxT )) || converged ) {
    // Gives access to velocity information on lattice
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    // Define output file
    std::string name =  "line.ghia.txt";
    //
    std::ofstream fout( name, std::ios::trunc );
    if ( !fout ) {
      clout << "Error: could not open " << name << std::endl;
    }
    // Interpolation functor with velocity information on lattice
    AnalyticalFfromSuperF3D<T> analytVel( velocity, true, 1 );
    Vector<int,17> y_coord( {128, 125, 124, 123, 122, 109, 94, 79, 64, 58, 36, 22, 13, 9, 8, 7, 0} );
    // ghia, ghia and shin, 1982: "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method";  Table 1, 2D data for comparison!
    Vector<T,17> vel_ghia_RE1000( {1.0,     0.65928, 0.57492, 0.51117, 0.46604,
                                   0.33304, 0.18719, 0.05702,-0.06080,-0.10648,
                                   -0.27805,-0.38289,-0.29730,-0.22220,-0.20196,
                                   -0.18109, 0.0
                                  } );
    Vector<T,17> vel_ghia_RE100( {1.0,     0.84123, 0.78871, 0.73722, 0.68717,
                                  0.23151, 0.00332,-0.13641,-0.20581,-0.21090,
                                  -0.15662,-0.10150,-0.06434,-0.04775,-0.04192,
                                  -0.03717, 0.0
                                 } );
    Vector<T,17> vel_simulation;
    // Define comparison values
    Vector<T,17> comparison = vel_ghia_RE1000;
    for ( int nY = 0; nY < 17; ++nY ) {
      // 17 data points evenly distributed between 0 and 1 (height)
      double postition[3] = {0.5, y_coord[nY]/128.0, 0.5};
      double velocity3[3] = {0.0};
      // Interpolate velocityField at "position"
      analytVel( velocity3, postition );
      // Save value of velocity (in x-direction) i
      vel_simulation[nY] = velocity3[0];
      // Write position and x-velocity in file
      fout << postition[1] << ", " << vel_simulation[nY] << std::endl;
    }
    // Write errors in file
    fout << "\nL2(line) absolute error: " << ( vel_simulation - comparison ).norm() / 17. << std::endl;
    fout << "L2(line) relative error: " << ( vel_simulation - comparison ).norm() / comparison.norm() << std::endl;
    fout.close();
    // Print errors on console
    clout << "absoluteErrorL2(line)=" << ( vel_simulation - comparison ).norm() / 17.
          << "; relativeErrorL2(line)=" << ( vel_simulation - comparison ).norm() / comparison.norm() << std::endl;
  }

}



int main( int argc, char **argv ) {

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  LBconverter<T> converter(
    ( int ) 3,                             // dim
    ( T )   1./30./N,                      // latticeL_
    ( T )   1e-1/M,                        // latticeU_
    ( T )   1./1000.,                      // charNu_
    ( T )   1.,                            // charL_ = 1
    ( T )   1.                             // charU_ = 1
  );
  converter.print();
  writeLogFile( converter, "cavity3d" );

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a unit cube by an indicator
  Vector<T,3> origin( T( 0 ) );
  Vector<T,3> extend( converter.getCharL() );
  IndicatorCuboid3D<T> cube( extend,origin );

  // Instantiation of a cuboid geometry with weights
  int noCuboids = singleton::mpi().getSize();
  if ( noCuboids < 7 ) {
    noCuboids = 7;
  }
  CuboidGeometry3D<T> cuboidGeometry( cube, converter.getLatticeL(), noCuboids );

  // Instantiation of a load balancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a super geometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, cube, superGeometry );


  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getOmega(), instances::getBulkMomenta<T,DESCRIPTOR>() );

  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T,DESCRIPTOR,ConstRhoBGKdynamics<T,DESCRIPTOR> >( sBoundaryCondition );

  prepareLattice( converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  util::ValueTracer<T> converge( converter.numTimeSteps(interval), epsilon );

  Timer<T> timer( converter.numTimeSteps( maxT ), converter.numNodes( 1 )*converter.numNodes( 1 )*converter.numNodes( 1 ) );
  timer.start();

  for ( int iT = 0; iT <= converter.numTimeSteps( maxT ); ++iT ) {

    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( sLattice, converter, superGeometry, iT, timer, converge.hasConverged() );
      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, superGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, superGeometry, iT, timer, converge.hasConverged() );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}
