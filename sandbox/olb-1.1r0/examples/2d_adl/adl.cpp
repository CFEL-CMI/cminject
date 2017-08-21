/*  Lattice Boltzmann code for simulating gas in aerodynamic lens 2D, 
 * written in C++, using the OpenLB
 * By Muhamed Amin, CMI-CFEL 
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif

#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;


typedef double T;
#define DESCRIPTOR D2Q9Descriptor


// Parameters for the simulation setup

//const T U = 3.0;
const int N = 50;        // resolution of the model
const int M = 1;        // time discretization refinement
//const T Re = 30.;       // Reynolds number
const T L = 0.01/N;     // latticeL
const T lx1   = 0.050;    // length of first part
const T ly1   = 0.010;   // height of first part
const T lx2   = 0.001;   // length of narrow part in meter
const T ly2   = 0.003;    // height of narrow part in meter
const T lx3   = 0.005;    // length of third part   
const T ly3   = 0.002;   // height of third part
const T maxPhysT = 2.; // max. simulation time in s, SI unit


// Stores geometry information in form of material numbers
void prepareGeometry( LBconverter<T> const& converter,
                      SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> first_part( lx1,ly1 );
  Vector<T,2> origin1(0.,0.);
  IndicatorCuboid2D<T> adl1( first_part, origin1 );

  Vector<T,2> narrow( lx2,ly1-ly2 );
  Vector<T,2> origin2(lx1/2,ly2);
  IndicatorCuboid2D<T> adl2( narrow, origin2 );

  Vector<T,2> third_part( lx3,ly1-ly3 );
  Vector<T,2> origin3(lx1-lx3,ly3);
  IndicatorCuboid2D<T> adl3( third_part, origin3 );
  


  IndicatorIdentity2D<T> ADL( adl1-adl3-adl2 );

  // material numbers from zero to 2 inside geometry defined by indicator
  superGeometry.rename( 0,2,ADL );
  superGeometry.rename( 2,1,1,1 );

  Vector<T,2> extendBC( 0.005,ly1 );
  Vector<T,2> originBC;
  // Set material number for inflow
  IndicatorCuboid2D<T> inflow( extendBC, originBC );
  superGeometry.rename( 2,3,1,inflow );
  // Set material number for outflow
  originBC[0] = lx1;
  extendBC[1] = ly1-ly3;
  IndicatorCuboid2D<T> outflow( extendBC, originBC );
  superGeometry.rename( 2,4,1,outflow );
  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( LBconverter<T> const& converter,
                     SuperLattice2D<T,DESCRIPTOR>& sLattice,
                     Dynamics<T,DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& bc,
                     SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>() );
  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );
  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T,DESCRIPTOR>() );
  // Material=3 -->bulk dynamics (inflow)
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );
  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );

  // Setting of the boundary conditions
  bc.addVelocityBoundary( superGeometry, 3, converter.getOmega() );
  bc.addPressureBoundary( superGeometry, 4, converter.getOmega() );

  // Initial conditions
//  AnalyticalConst2D<T,T> ux( 5. );
//  AnalyticalConst2D<T,T> uy( 0. );
//  AnalyticalConst2D<T,T> rho( 0.002 );
//  AnalyticalComposed2D<T,T> u( ux,uy );

  //Initialize all values of distribution functions to their local equilibrium
//  sLattice.defineRhoU( superGeometry, 1, rho, u );
//  sLattice.iniEquilibrium( superGeometry, 1, rho, u );
//    sLattice.defineRhoU( superGeometry, 3, rho, u );
//    sLattice.iniEquilibrium( superGeometry, 3, rho, u );
//  sLattice.defineRhoU( superGeometry, 4, rho, u );
//  sLattice.iniEquilibrium( superGeometry, 4, rho, u );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( LBconverter<T> const& converter,
                        SuperLattice2D<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // time for smooth start-up
  int iTmaxStart = converter.numTimeSteps( maxPhysT*0.2 );
  int iTupdate = 5;

  if ( iT%iTupdate == 0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));
    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );
    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getLatticeU()*3./2.*frac[0];
    T distance2Wall = converter.getLatticeL()/2.;
 
//    cout<<"maxVelocity   :"<<maxVelocity<<endl; 

    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );
    // define lattice speed on inflow
    sLattice.defineU( superGeometry, 3, poiseuilleU );
  }
}

// write data to termimal and file system
void getResults( SuperLattice2D<T,DESCRIPTOR>& sLattice,
                 LBconverter<T> const& converter, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer ) {
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "bstep2d" );

  if ( iT==0 ) {
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeGeometry2D<T,DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T,DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T,DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes every 0.2  seconds
  if ( iT%converter.numTimeSteps( 0.2 )==0 ) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    // write vtk to file system
    vtmWriter.write( iT );
    sLattice.communicate();
    SuperEuklidNorm2D<T,DESCRIPTOR> normVel( velocity );
    BlockLatticeReduction2D<T,DESCRIPTOR> planeReduction( normVel );
    BlockGifWriter<T> gifWriter;
    // write image to file system
    gifWriter.write( planeReduction, iT, "vel" );


    AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> intpolatePressure( pressure, true );
    AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> intpolateVelocity( velocity, true );
    T point1[2] = {};
    point1[0] = 0.001;
    point1[1] = 0.005;
    T p1;
    T v1;
    intpolatePressure( &p1,point1 );
    intpolateVelocity( &v1,point1 );
    clout << "######################  The pressure1 at ( "<<point1[0]<<", "<<point1[1]<<") is: " << p1 <<" ########################"<<endl;
    clout << "######################  The velocity1 at ( "<<point1[0]<<", "<<point1[1]<<") is: " << v1 <<" ########################"<<endl;
  }

  // Writes every 0.1 simulated
  if ( iT%converter.numTimeSteps( 0.1 )==0 ) {
    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.physTime( iT ) );
  }

}


int main( int argc, char* argv[] ) {
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );  // set output directory
  OstreamManager clout( std::cout, "main" );

  LBconverter<T> converter(
    ( int ) 2,                             // dim
    ( T )   L,                             // latticeL_
    ( T )   0.01/M,                        // latticeU_
    ( T )   0.0005,                  // charNu_
    ( T )   0.01,                         // charL_ = 1
    ( T )   2.0,                            // charU_ = 1
    ( T )   1., 
    ( T )   0. 
  );
  converter.print();

  writeLogFile( converter, "bstep2d" );

  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( lx1, ly1);
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getLatticeL(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );


  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T,DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T,DESCRIPTOR> bulkDynamics (
    converter.getOmega(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  // createInterpBoundaryCondition2D<T,DESCRIPTOR>(sBoundaryCondition);
  createLocalBoundaryCondition2D<T,DESCRIPTOR>( sBoundaryCondition );

  prepareLattice( converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  Timer<T> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( int iT = 0; iT < converter.numTimeSteps( maxPhysT ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}