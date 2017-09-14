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


/*
void write_v_and_p(SuperLatticePhysVelocity2D<T,DESCRIPTOR>& velocity,
                   SuperLatticePhysPressure2D<T,DESCRIPTOR>& pressure,
                   T x, T y, T res ) {

    AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> intpolatePressure( pressure, true );
    AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> intpolateVelocity( velocity, true );
    T point[3] = {};
    point[0] = 0.000;
    point[1] = 0.000;
    point[2]=  0.000;
    T p;
    T v;

    ofstream myfile;
    myfile.open ("data.txt");
    myfile << "X Y V P "<<endl;
    for (int i=0; i<(x/res); i++) {
     for (int j=0; j<(y/res); j++) {
         point[1] = point[1] + res;
         intpolatePressure( &p,point );
         intpolateVelocity( &v,point );
         myfile << point[0] <<" "<< point[1] << " " << p <<" "<< v<<endl;;
     }
         point[0] = point[0] + res;
         point[1]= 0;
         cout<<"#";
    }
   myfile.close();


}
*/



const int N = 100;    // resolution of the model
const int M = 1;    // time discretization refinement
T maxPhysT = 200.0; // max. simulation time in s, SI unit

void prepareGeometry( LBconverter<T> const& converter, IndicatorF3D<T>& indicator, SuperGeometry3D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;


  superGeometry.rename( 0,2,indicator );

  // Definition of the geometry of the venturi
  Vector<T,3> C0( 0,0,0 );
  Vector<T,3> C1( 0.01,0,0 );
  Vector<T,3> C2( 0.1,0,0 );
  Vector<T,3> C3( 0.12,0,0 );
  Vector<T,3> C4( 0.2,0,0 );
  Vector<T,3> C5( 0.24,0,0 );
  Vector<T,3> C6( 0.25,0,0 );

  T radius1 = 0.01;  // radius of the inlet 
  T radius2 = 0.0025;  // radius of the lens 
  T radius3 = 0.001;   // radius of the small exit

  IndicatorCylinder3D<T> inflow( C0, C1, radius1 );
  IndicatorCylinder3D<T> cyl1( C1, C2, radius1 );
  IndicatorCylinder3D<T> cyl2( C2, C3, radius2 );
  IndicatorCylinder3D<T> cyl3( C3, C4, radius1 );
  IndicatorCylinder3D<T> cyl4( C4, C5, radius3 );
  IndicatorCylinder3D<T> outflow( C5, C6, radius3 );



  // Set boundary voxels by rename material numbers

  superGeometry.rename( 2,1,cyl1 );
  superGeometry.rename( 2,1,cyl2 );
  superGeometry.rename( 2,1,cyl3 );
  superGeometry.rename( 2,1,cyl4 );
  superGeometry.rename( 1,3,inflow );
  superGeometry.rename( 1,4,outflow );


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
                     sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
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

  // Initialize all values of distribution functions to their local equilibrium
  // Initial conditions

  // Lattice initialize
  sLattice.initialize();



  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing sinuidal inflow for the first iTMax timesteps
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        LBconverter<T> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.numTimeSteps( maxPhysT*0.8 );
  int iTperiod = 50;

  if ( iT==0 ) {
    // Make the lattice ready for simulation
    sLattice.initialize();
  }

  else if ( iT%iTperiod==0 && iT<= iTmaxStart ) {
    //clout << "Set Boundary Values ..." << std::endl;

    //SinusStartScale<T,int> startScale(iTmaxStart, (T) 1);
    PolynomialStartScale<T,int> startScale( iTmaxStart, T( 1 ) );
    int iTvec[1]= {iT};
    T frac = T();
    startScale( &frac,iTvec );

    // Creates and sets the Poiseuille inflow profile using functors
    CirclePoiseuille3D<T> poiseuilleU( superGeometry, 3, frac*converter.getLatticeU(), converter.getLatticeL() );
    sLattice.defineU( superGeometry, 3, poiseuilleU );

  //  clout << "step=" << iT << "; scalingFactor=" << frac <<"; lattice spead "<<converter.getLatticeU()<< std::endl;
  }
 // clout << "Set Boundary Values ... ok" << std::endl;
}

void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 LBconverter<T>& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "venturi3d" );

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
  if ( iT%converter.numTimeSteps( 1. )==0 ) {
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
  }

  // Writes output on the console
  if ( iT%converter.numTimeSteps( 1. )==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT, converter.physTime( iT ) );

  }
}



int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  // clout.setMultiOutput(true);


  LBconverter<T> converter(
    ( int ) 3,                             // dim
    ( T )   0.1/N,                          // latticeL_
    ( T )   0.01,                        // latticeU_
    ( T )   0.12,                           // charNu_
    ( T )   0.01,                           // charL_ = 1
    ( T )   10.,                             // charU_ = 1
    ( T )   0.000164
  );


  converter.print();
  writeLogFile( converter, "venturi3d" );

  // === 2nd Step: Prepare Geometry ===


  Vector<T,3> origin;
  Vector<T,3> extend( 0.25, 0.005+2.*converter.getLatticeL(), 0.005+2.*converter.getLatticeL() );

  IndicatorCuboid3D<T> cuboid( extend,origin );

  CuboidGeometry3D<T> cuboidGeometry( cuboid, converter.getLatticeL(), singleton::mpi().getSize() );
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // === 2nd Step: Prepare Geometry ===

 
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );
  prepareGeometry( converter, cuboid, superGeometry );


  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  RLBdynamics<T, DESCRIPTOR> bulkDynamics( converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>() );

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T, DESCRIPTOR> ( sBoundaryCondition );

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );

  prepareLattice( sLattice, converter, bulkDynamics, sBoundaryCondition, sOffBoundaryCondition, superGeometry );

  Timer<T> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  getResults( sLattice, converter, 0, superGeometry, timer );

  // === 4th Step: Main Loop with Timer ===
  for ( int iT = 0; iT <= converter.numTimeSteps( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}



