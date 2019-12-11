/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014 Mathias J. Krause, Thomas Henn,
 *  Cyril Masquelier
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

/* venturi3d.cpp:
 * This example examines a steady flow in a venturi tube. At the
 * main inlet, a Poiseuille profile is imposed as Dirichlet velocity
 * boundary condition, whereas at the outlet and the minor inlet
 * a Dirichlet pressure condition is set by p=0 (i.e. rho=1).
 *
 * The example shows the usage of the Indicator functors to
 * build up a geometry and explains how to set boundary conditions
 * automatically.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
#include "olb3D.hh"     // Include full template code
#endif

#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor





//const int N = 100;    // resolution of the model
//const int M = 1;    // time discretization refinement
T maxPhysT = 150.0; // max. simulation time in s, SI unit

void prepareGeometry( LBconverter<T> const& converter, IndicatorF3D<T>& indicator, SuperGeometry3D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;


  superGeometry.rename( 0,2,indicator );
  const T L = 0.015 + converter.getLatticeL();
  const T e = converter.getLatticeL();

  Vector<T,3> Ci0( 0.000,L,L + L/2 );
  Vector<T,3> Ci1( 0.001,L,L + L/2 );
  Vector<T,3> Cie( e,L,L + L/2 );

  Vector<T,3> C0( 0.000,L,L );
  Vector<T,3> Co0( e,L,L );
  Vector<T,3> C1( 0.001,L,L );
  Vector<T,3> C2( 0.002,L,L );
  Vector<T,3> C3( 0.003,L,L );
  Vector<T,3> C4( 0.013,L,L );
  Vector<T,3> C5( 0.033,L,L );
  Vector<T,3> C6( 0.043,L,L );
  Vector<T,3> C7( 0.044,L,L );
  Vector<T,3> Co1( 0.044-e,L,L );

  T radius1 = 0.001;  // radius of the inlet 
  T radius2 = 0.015;  // radius of the lens 
  T radius3 = 0.002;   // radius of the small exit


  IndicatorCylinder3D<T> inflow( Ci0, Cie, radius1 );

  IndicatorCylinder3D<T> cyl1( Cie, Ci1, radius1 );
  IndicatorCylinder3D<T> cyl2( C0, C1, radius1 );
  IndicatorCylinder3D<T> cyl3( C1, C2, radius2 );
  IndicatorCylinder3D<T> cyl4( C2, C3, radius3 );
  IndicatorCone3D<T> co1( C3, C4, radius3, radius2 );
  IndicatorCylinder3D<T> cyl5( C4, C5, radius2 );
  IndicatorCone3D<T> co2( C5, C6, radius2, radius1 );
  IndicatorCylinder3D<T> cyl6( C6, C7, radius1 );

  IndicatorCylinder3D<T> outflow0( C0, Co0, radius3 );
  IndicatorCylinder3D<T> outflow1( Co1, C7, radius3 );


  // Set boundary voxels by rename material numbers

  superGeometry.rename( 2,1,cyl1 );
  superGeometry.rename( 2,1,cyl2 );
  superGeometry.rename( 2,1,cyl3 );
  superGeometry.rename( 2,1,cyl4 );
  superGeometry.rename( 2,1,cyl5 );
  superGeometry.rename( 2,1,cyl6 );
  superGeometry.rename( 2,1,co1 );
  superGeometry.rename( 2,1,co2 );
  superGeometry.rename( 2,3,inflow );
  superGeometry.rename( 2,4,outflow0 );
  superGeometry.rename( 2,4,outflow1 );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();


  superGeometry.getStatistics().print();
  superGeometry.communicate();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     LBconverter<T> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     SuperGeometry3D<T>& superGeometry ) {



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
  bc.addVelocityBoundary( superGeometry, 3, omega );
  bc.addPressureBoundary( superGeometry, 4, omega );

  sLattice.initialize();



  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing sinuidal inflow for the first iTMax timesteps
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        LBconverter<T> const& converter,
                        SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );
    std::vector<T> velocity( 3,T( 0 ) );
    velocity[0] = converter.getLatticeU();
    AnalyticalConst3D<T,T> uF( velocity );
    sLattice.defineU( superGeometry, 3, uF );  
}

void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 LBconverter<T>& converter, int iT, bool converged,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer, float convEach ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "bgc" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the vtm files
  if ( converged ) {
    // Create the data-reading functors...
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockLatticeReduction3D<T, DESCRIPTOR> planeReduction( normVel, 0, 0, -1 );
    BlockGifWriter<T> gifWriter;
    gifWriter.write( planeReduction, iT, "vel" ); // scaled
    //write_v_and_p(velocity, pressure, 0.05, 0.02, 0.00025, converter);

  }

  // Writes output on the console
  if ( iT%converter.numTimeSteps( convEach )==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT, converter.physTime( iT ) );

  }
}



int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  // clout.setMultiOutput(true);

  clout<<argv[1]<<endl;
  float InletSpeed = atof(argv[1]);
  float KinematicV = atof(argv[2]);
  float Density = atof(argv[3]);
  double Conv = atof(argv[4]);
  float U = atof(argv[5]);
  std::string path = argv[6];
  float convEach = U/(75*KinematicV*50);

  LBconverter<T> converter(
    ( int ) 3,                             // dim
    ( T )   0.0004,                        // latticeL_
    ( T )   U,                             // latticeU_
    ( T )   KinematicV,                      // charNui_ for He at 4K =0.0015
    ( T )   0.004,                         // charL_ = 1
    ( T )   InletSpeed,                    // charU_ = 1 
    ( T )   Density                        // 0.0164 for He at 4K
  );

  singleton::directories().setOutputDir( path+"/LBM/" );
  converter.print();
  writeLogFile( converter, "bgc" );

  // === 2nd Step: Prepare Geometry ===

  Vector<T,3> origin;
  Vector<T,3> extend( 0.045+converter.getLatticeL(), 0.03+2.*converter.getLatticeL(), 0.03+2.*converter.getLatticeL() );

  IndicatorCuboid3D<T> cuboid( extend,origin );

  CuboidGeometry3D<T>* cuboidGeometry = new CuboidGeometry3D<T>( cuboid, converter.getLatticeL(), singleton::mpi().getSize());

  //CuboidGeometry3D<T> cuboidGeometry( cuboid, converter.getLatticeL(), singleton::mpi().getSize() );

  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>( *cuboidGeometry );

  //HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // === 2nd Step: Prepare Geometry ===

 
  SuperGeometry3D<T> superGeometry( *cuboidGeometry, *loadBalancer, 2 );
  prepareGeometry( converter, cuboid, superGeometry );


  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>() );

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T, DESCRIPTOR> ( sBoundaryCondition );

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );

  prepareLattice( sLattice, converter, bulkDynamics, sBoundaryCondition, superGeometry );

  Timer<T> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  getResults( sLattice, converter, 0, false, superGeometry, timer, convEach );
  setBoundaryValues( sLattice, converter, superGeometry );

  util::ValueTracer<T> converge( converter.numTimeSteps(convEach), Conv );
  // === 4th Step: Main Loop with Timer ===
  for ( int iT = 0; iT <= converter.numTimeSteps( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, false, superGeometry, timer, convEach );

    // === check for convergence 
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true);
    if (converge.hasConverged()) {
    clout << " Simulation converged ... Writing the Flow Field " << endl; 
    getResults( sLattice, converter, iT, true, superGeometry, timer, convEach );
    break;
    }
  }

  timer.stop();
  timer.printSummary();
}

