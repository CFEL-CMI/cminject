/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2014 Mathias J. Krause
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

/* aorta3d.cpp:
 * In this example the fluid flow through a bifurcation is
 * simulated. The geometry is obtained from a mesh in stl-format.
 * With Bouzidi boundary conditions the curved boundary is
 * adequately mapped and initialized fully automatically. As
 * dynamics a Smagorinsky turbulent BGK model is used to stabilize
 * the simulation for low resolutions. As output the flux at the
 * inflow and outflow region is computed. The results has been
 * validated by comparison with other results obtained with FEM
 * and FVM.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code;
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor


//simulation parameters
const int N = 1;             // resolution of the model
const int M = 1;             // time discretization refinement
const bool bouzidiOn = true; // choice of boundary condition
const T maxPhysT = 2.;       // max. simulation time in s, SI unit


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( LBconverter<T> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,stlReader );

  superGeometry.clean();

  // Set material number for inflow
  IndicatorCircle3D<T> inflow(  0.218125 ,0.249987 ,0.0234818, 0., 1.,0., 0.0112342 );
  IndicatorCylinder3D<T> layerInflow( inflow, 2.*converter.getLatticeL() );
  superGeometry.rename( 2,3,1,layerInflow );

  // Set material number for outflow0
  //IndicatorCircle3D<T> outflow0(0.2053696,0.0900099,0.0346537,  2.5522,5.0294,-1.5237, 0.0054686 );
  IndicatorCircle3D<T> outflow0( 0.2053696,0.0900099,0.0346537, 0.,-1.,0., 0.0054686 );
  IndicatorCylinder3D<T> layerOutflow0( outflow0, 2.*converter.getLatticeL() );
  superGeometry.rename( 2,4,1,layerOutflow0 );

  // Set material number for outflow1
  //IndicatorCircle3D<T> outflow1(0.2388403,0.0900099,0.0343228, -1.5129,5.1039,-2.8431, 0.0058006 );
  IndicatorCircle3D<T> outflow1( 0.2388403,0.0900099,0.0343228, 0.,-1.,0., 0.0058006 );
  IndicatorCylinder3D<T> layerOutflow1( outflow1, 2.*converter.getLatticeL() );
  superGeometry.rename( 2,5,1,layerOutflow1 );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean( 3 );
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& lattice,
                     LBconverter<T> const& converter, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& bc,
                     sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                     STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getOmega();

  // material=0 --> do nothing
  lattice.defineDynamics( superGeometry,0,&instances::getNoDynamics<T, DESCRIPTOR>() );

  // material=1 --> bulk dynamics
  lattice.defineDynamics( superGeometry,1,&bulkDynamics );

  if ( bouzidiOn ) {
    // material=2 --> no dynamics + bouzidi zero velocity
    lattice.defineDynamics( superGeometry,2,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addZeroVelocityBoundary( superGeometry,2,stlReader );
    // material=3 --> no dynamics + bouzidi velocity (inflow)
    lattice.defineDynamics( superGeometry,3,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addVelocityBoundary( superGeometry,3,stlReader );
  } else {
    // material=2 --> bounceBack dynamics
    lattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );
    // material=3 --> bulk dynamics + velocity (inflow)
    lattice.defineDynamics( superGeometry,3,&bulkDynamics );
    bc.addVelocityBoundary( superGeometry,3,omega );
  }

  // material=4,5 --> bulk dynamics + pressure (outflow)
  lattice.defineDynamics( superGeometry,4,&bulkDynamics );
  lattice.defineDynamics( superGeometry,5,&bulkDynamics );
  bc.addPressureBoundary( superGeometry,4,omega );
  bc.addPressureBoundary( superGeometry,5,omega );

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  lattice.defineRhoU( superGeometry,1,rhoF,uF );
  lattice.iniEquilibrium( superGeometry,1,rhoF,uF );
  lattice.defineRhoU( superGeometry,3,rhoF,uF );
  lattice.iniEquilibrium( superGeometry,3,rhoF,uF );
  lattice.defineRhoU( superGeometry,4,rhoF,uF );
  lattice.iniEquilibrium( superGeometry,4,rhoF,uF );
  lattice.defineRhoU( superGeometry,5,rhoF,uF );
  lattice.iniEquilibrium( superGeometry,5,rhoF,uF );

  // Lattice initialize
  lattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing sinuidal inflow
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                        LBconverter<T> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry ) {

  // No of time steps for smooth start-up
  int iTperiod = converter.numTimeSteps( 0.5 );
  int iTupdate = 50;

  if ( iT%iTupdate == 0 ) {
    // Smooth start curve, sinus
    SinusStartScale<T,int> nSinusStartScale( iTperiod,converter.getLatticeU() );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1]= {iT};
    T maxVelocity[1]= {T()};
    nSinusStartScale( maxVelocity,iTvec );
    CirclePoiseuille3D<T> velocity( superGeometry,3,maxVelocity[0] );

    if ( bouzidiOn ) {
      offBc.defineU( superGeometry,3,velocity );
    } else {
      sLattice.defineU( superGeometry,3,velocity );
    }
  }
}

// Computes flux at inflow and outflow
void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 LBconverter<T>& converter, int iT,
                 Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer, STLreader<T>& stlReader ) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "aorta3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int vtkIter  = converter.numTimeSteps( .1 );
  const int statIter = converter.numTimeSteps( .1 );

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

  // Writes the vtk files
  if ( iT%vtkIter==0 ) {
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockLatticeReduction3D<T, DESCRIPTOR> planeReduction( normVel, 0, 0, -1 );
    BlockGifWriter<T> gifWriter;
    //gifWriter.write(planeReduction, 0, 0.7, iT, "vel"); //static scale
    gifWriter.write( planeReduction, iT, "vel" ); // scaled
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.physTime( iT ) );

    // Flux at the inflow and outflow region
    std::list<int> materials;
    materials.push_back( 1 );
    materials.push_back( 3 );
    materials.push_back( 4 );
    materials.push_back( 5 );

    IndicatorCircle3D<T> inflow(  0.218125 ,0.249987-2.*converter.getLatticeL() ,0.0234818, 0., -1.,0., 0.0112342+2*converter.getLatticeL() );
    SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR> vFluxInflow( sLattice, converter, superGeometry, inflow, materials );
    vFluxInflow.print( "inflow","ml/s" );
    SuperLatticePhysPressureFlux3D<T, DESCRIPTOR> pFluxInflow( sLattice, converter, superGeometry, inflow, materials );
    pFluxInflow.print( "inflow","N","mmHg" );

    IndicatorCircle3D<T> outflow0( 0.2053696,0.0900099+2.*converter.getLatticeL(),0.0346537, 0.,1.,0., 0.0054686+2*converter.getLatticeL() );
    SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR> vFluxOutflow0( sLattice, converter, superGeometry, outflow0, materials );
    vFluxOutflow0.print( "outflow0","ml/s" );
    SuperLatticePhysPressureFlux3D<T, DESCRIPTOR> pFluxOutflow0( sLattice, converter, superGeometry, outflow0, materials );
    pFluxOutflow0.print( "outflow0","N","mmHg" );

    IndicatorCircle3D<T> outflow1( 0.2388403,0.0900099+2.*converter.getLatticeL(),0.0343228, 0.,1.,0., 0.0058006+2*converter.getLatticeL() );
    SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR> vFluxOutflow1( sLattice, converter, superGeometry, outflow1, materials );
    vFluxOutflow1.print( "outflow1","ml/s" );
    SuperLatticePhysPressureFlux3D<T, DESCRIPTOR> pFluxOutflow1( sLattice, converter, superGeometry, outflow1, materials );
    pFluxOutflow1.print( "outflow1","N","mmHg" );

    if ( bouzidiOn ) {
      SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, superGeometry, stlReader, 3 );
      SuperMax3D<T> yPlusMaxF( yPlus, superGeometry, 1 );
      int input[4]= {};
      T yPlusMax = T();
      yPlusMaxF( &yPlusMax,input );
      clout << "yPlusMax=" << yPlusMax << std::endl;
    }
  }

  if ( sLattice.getStatistics().getMaxU() > 0.3 ) {
    clout << "PROBLEM uMax=" << sLattice.getStatistics().getMaxU() << std::endl;
    vtmWriter.write( iT );
    std::exit( 0 );
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
    3,                                     // dim
    ( T ) 0.00056179/N,                    // latticeL_
    ( T ) 0.05/M,                          // latticeU_
    ( T ) 0.003/1055.,                     // charNu_
    ( T ) 0.02246,                         // charL_ = 1,
    ( T ) 0.45,                            // charU_ = 1,
    ( T ) 1055,                            // charRho_ = 1,
    ( T ) 0                                // charPressure_ = 0
  );
  converter.print();
  writeLogFile( converter, "aorta3d" );


  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  STLreader<T> stlReader( "aorta3d.stl", converter.getLatticeL(), 0.001, 0, true );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getLatticeL() );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = std::min( 16*N,2*singleton::mpi().getSize() );
#else
  const int noOfCuboids = 2;
#endif
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getLatticeL(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  SmagorinskyBGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getOmega(),
      instances::getBulkMomenta<T, DESCRIPTOR>(), 0.1, converter.getLatticeL(), converter.physTime() );
  //BGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getOmega(),
  //                                        instances::getBulkMomenta<T, DESCRIPTOR>());

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryCondition );
  // createLocalBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );

  Timer<T> timer1( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer1.start();

  prepareLattice( sLattice, converter, bulkDynamics,
                  sBoundaryCondition, sOffBoundaryCondition,
                  stlReader, superGeometry );

  timer1.stop();
  timer1.printSummary();

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  Timer<T> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( int iT = 0; iT <= converter.numTimeSteps( maxPhysT ); iT++ ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, sOffBoundaryCondition, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, bulkDynamics, superGeometry, timer, stlReader );
  }

  timer.stop();
  timer.printSummary();
}
