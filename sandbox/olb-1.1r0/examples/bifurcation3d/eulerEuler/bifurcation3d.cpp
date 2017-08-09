/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2016 Robin Trunk, Mathias J. Krause
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

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation by imposing
 * a inflow condition with rho = 1 on the second euler phase at the inlet.
 * The particles are simulated as a continuum with a advection-diffusion equation
 * and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5
 *  *
 */

#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!

using namespace std;
using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;
#define NSDESCRIPTOR D3Q19Descriptor
#define ADDESCRIPTOR particleAdvectionDiffusionD3Q7Descriptor

const T Re = 50;               // Reynolds number
const int N = 1;               // resolution of the model
const int M = 1;               // time discretization refinement
const int iTperiod = 100;      // amount of timesteps when new boundary conditions are reset and results are visualized
const T diffusion = 1.e-6;     // diffusion coefficient for advection-diffusion equation
const T radius = 1.5e-4;       // particles radius
const T partRho = 998.2;       // particles density
const T maxPhysT = 10.;        // max. simulation time in s, SI unit

// center of inflow and outflow regions [m]
Vector<T,3> inletCenter( T(), T(), 0.0786395 );
Vector<T,3> outletCenter0( -0.0235929682287551, -0.000052820468762797, -0.021445708949909 );
Vector<T,3> outletCenter1( 0.0233643529416147, 0.00000212439067050152, -0.0211994104877918 );

// radii of inflow and outflow regions [m]
T inletRadius = 0.00999839;
T outletRadius0 = 0.007927;
T outletRadius1 = 0.00787134;

// normals of inflow and outflow regions
Vector<T,3> inletNormal( T(), T(), T( -1 ) );
Vector<T,3> outletNormal0( 0.505126, -0.04177, 0.862034 );
Vector<T,3> outletNormal1( -0.483331, -0.0102764, 0.875377 );

void prepareGeometry( LBconverter<T> const& converter,
                      IndicatorF3D<T>& indicator, STLreader<T>& stlReader,
                      SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout, "prepareGeometry" );

  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 2, indicator );
  superGeometry.rename( 2, 1, stlReader );

  superGeometry.clean();

  // rename the material at the inlet
  IndicatorCircle3D<T> inletCircle( inletCenter, inletNormal, inletRadius );
  IndicatorCylinder3D<T> inlet( inletCircle, 2 * converter.getLatticeL() );
  superGeometry.rename( 2, 3, 1, inlet );

  // rename the material at the outlet0
  IndicatorCircle3D<T> outletCircle0( outletCenter0, outletNormal0, 0.95*outletRadius0 );
  IndicatorCylinder3D<T> outlet0( outletCircle0, 4 * converter.getLatticeL() );
  superGeometry.rename( 2, 4, outlet0 );

  // rename the material at the outlet1
  IndicatorCircle3D<T> outletCircle1( outletCenter1, outletNormal1, 0.95*outletRadius1 );
  IndicatorCylinder3D<T> outlet1( outletCircle1, 4 * converter.getLatticeL() );
  superGeometry.rename( 2, 5, outlet1 );

  superGeometry.clean();
  superGeometry.innerClean( true );
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice( SuperLattice3D<T, NSDESCRIPTOR>& sLatticeNS,
                     SuperLattice3D<T, ADDESCRIPTOR>& sLatticeAD,
                     LBconverter<T> const& converter,
                     Dynamics<T, NSDESCRIPTOR>& bulkDynamics,
                     Dynamics<T, ADDESCRIPTOR>& bulkDynamicsAD,
                     sOnLatticeBoundaryCondition3D<T, NSDESCRIPTOR>& bc,
                     sOnLatticeBoundaryCondition3D<T, ADDESCRIPTOR>& bcAD,
                     SuperGeometry3D<T>& superGeometry,
                     T omegaAD ) {

  OstreamManager clout( std::cout, "prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getOmega();

  // Material=0 --> do nothing
  sLatticeNS.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>() );
  sLatticeAD.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, ADDESCRIPTOR>() );

  // Material=1 --> bulk dynamics
  sLatticeNS.defineDynamics( superGeometry, 1, &bulkDynamics );
  sLatticeAD.defineDynamics( superGeometry, 1, &bulkDynamicsAD );

  // Material=2 --> bounce-back /
  sLatticeNS.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, NSDESCRIPTOR>() );
  sLatticeAD.defineDynamics( superGeometry, 2, &instances::getZeroDistributionDynamics<T, ADDESCRIPTOR>() );

  // Material=3 --> bulk dynamics (inflow)
  sLatticeNS.defineDynamics( superGeometry, 3, &bulkDynamics );
  sLatticeAD.defineDynamics( superGeometry, 3, &bulkDynamicsAD );

  // Material=4,5 -->bulk dynamics / do-nothing (outflow)
  sLatticeNS.defineDynamics( superGeometry, 4, &bulkDynamics );
  sLatticeAD.defineDynamics( superGeometry, 4, &instances::getNoDynamics<T, ADDESCRIPTOR>() );
  sLatticeNS.defineDynamics( superGeometry, 5, &bulkDynamics );
  sLatticeAD.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T, ADDESCRIPTOR>() );

  // Setting of the boundary conditions
  bc.addPressureBoundary( superGeometry, 3, omega );
  bc.addVelocityBoundary( superGeometry, 4, omega );
  bc.addVelocityBoundary( superGeometry, 5, omega );
  bcAD.addZeroDistributionBoundary( superGeometry, 2 );
  bcAD.addTemperatureBoundary( superGeometry, 3, omegaAD );
  bcAD.addConvectionBoundary( superGeometry, 4 );
  bcAD.addConvectionBoundary( superGeometry, 5 );
  bcAD.addExtFieldBoundary( superGeometry, 2, ADDESCRIPTOR<T>::ExternalField::velocityBeginsAt );
  bcAD.addExtFieldBoundary( superGeometry, 3, ADDESCRIPTOR<T>::ExternalField::velocityBeginsAt );
  bcAD.addExtFieldBoundary( superGeometry, 4, ADDESCRIPTOR<T>::ExternalField::velocityBeginsAt );
  bcAD.addExtFieldBoundary( superGeometry, 5, ADDESCRIPTOR<T>::ExternalField::velocityBeginsAt );

  // Initial conditions
  AnalyticalConst3D<T,T> rho1( 1. );
  AnalyticalConst3D<T,T> rho0( 1.e-8 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> u0( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLatticeNS.defineRhoU( superGeometry,1,rho1,u0 );
  sLatticeNS.iniEquilibrium( superGeometry,1,rho1,u0 );
  sLatticeNS.defineRhoU( superGeometry,2,rho1,u0 );
  sLatticeNS.iniEquilibrium( superGeometry,2,rho1,u0 );
  sLatticeNS.defineRhoU( superGeometry,3,rho1,u0 );
  sLatticeNS.iniEquilibrium( superGeometry,3,rho1,u0 );
  sLatticeNS.defineRhoU( superGeometry,4,rho1,u0 );
  sLatticeNS.iniEquilibrium( superGeometry,4,rho1,u0 );
  sLatticeNS.defineRhoU( superGeometry,5,rho1,u0 );
  sLatticeNS.iniEquilibrium( superGeometry,5,rho1,u0 );
  sLatticeAD.defineRho( superGeometry,3,rho1 );
  sLatticeAD.iniEquilibrium( superGeometry,1,rho0,u0 );
  sLatticeAD.iniEquilibrium( superGeometry,2,rho0,u0 );
  sLatticeAD.iniEquilibrium( superGeometry,4,rho0,u0 );
  sLatticeAD.iniEquilibrium( superGeometry,5,rho0,u0 );

  // Lattice initialize
  sLatticeNS.initialize();
  sLatticeAD.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

void setBoundaryValues( SuperLattice3D<T, NSDESCRIPTOR>& sLatticeNS,
                        LBconverter<T> const& converter, int iT, T maxPhysT,
                        SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout, "setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.numTimeSteps( 0.8*maxPhysT );
  // Set inflow velocity
  T maxVelocity = converter.getLatticeU() * 3. / 4. * std::pow( inletRadius, 2 ) / std::pow( outletRadius0, 2 );
  if ( iT % iTperiod == 0 ) {
    if ( iT <= iTmaxStart ) {
      SinusStartScale<T, int> startScale( iTmaxStart, T( 1 ) );
      int iTvec[1] = { iT };
      T frac[1] = { T( 0 ) };
      startScale( frac, iTvec );
      maxVelocity *= frac[0];
    }

    CirclePoiseuille3D<T> poiseuilleU4( outletCenter0[0], outletCenter0[1],
                                        outletCenter0[2], outletNormal0[0],
                                        outletNormal0[1], outletNormal0[2],
                                        outletRadius0 * 0.95, -maxVelocity );

    CirclePoiseuille3D<T> poiseuilleU5( outletCenter1[0], outletCenter1[1],
                                        outletCenter1[2], outletNormal1[0],
                                        outletNormal1[1], outletNormal1[2],
                                        outletRadius1 * 0.95, -maxVelocity );

    sLatticeNS.defineU( superGeometry, 4, poiseuilleU4 );
    sLatticeNS.defineU( superGeometry, 5, poiseuilleU5 );
  }
}

void getResults( SuperLattice3D<T, NSDESCRIPTOR>& sLatticeNS,
                 SuperLattice3D<T, ADDESCRIPTOR>& sLatticeAD,
                 LBconverter<T>& converter, int iT, SuperGeometry3D<T>& superGeometry,
                 Timer<double>& timer ) {

  OstreamManager clout( std::cout, "getResults" );
  SuperVTMwriter3D<T> vtmWriter( "bifurcation3d_fluid" );
  SuperVTMwriter3D<T> vtmWriterAD( "bifurcation3d_particle" );
  SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );
  SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure( sLatticeNS, converter );
  SuperLatticeDensity3D<T, ADDESCRIPTOR> particles( sLatticeAD );
  SuperLatticePhysExternal3D<T, ADDESCRIPTOR> extField( sLatticeAD, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriterAD.addFunctor( particles );
  vtmWriterAD.addFunctor( extField );

  if ( iT == 0 ) {
    writeLogFile( converter, "bifurcation3d" );
    SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry( sLatticeNS, superGeometry );
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid( sLatticeNS );
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank( sLatticeNS );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
    vtmWriterAD.createMasterFile();

    // Print some output of the chosen simulation setup
    clout << "N=" << N << "; maxTimeSteps=" << converter.numTimeSteps( maxPhysT )
          << "; noOfCuboid=" << superGeometry.getCuboidGeometry().getNc() << "; Re=" << Re
          << "; St=" << ( 2.*partRho*radius*radius*converter.getCharU() ) / ( 9.*converter.getCharNu()*converter.getCharRho()*converter.getCharL() )
          << std::endl;
  }

  if ( iT % iTperiod == 0 ) {
    // Writes the vtk files
    vtmWriter.write( iT );
    vtmWriterAD.write( iT );

    // GIF Writer
    SuperEuklidNorm3D<T, NSDESCRIPTOR> normVel( velocity );
    BlockLatticeReduction3D<T, NSDESCRIPTOR> planeReductionVelocity( normVel, 0, 1, 0 );
    BlockLatticeReduction3D<T, ADDESCRIPTOR> planeReductionParticles( particles, 0, 1, 0 );
    BlockGifWriter<T> gifWriter;
    gifWriter.write( planeReductionVelocity, 0, 0.06, iT, "vel" ); //static scale
    gifWriter.write( planeReductionParticles, iT, "particles" ); // scaled

    // Writes output on the console
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print( iT, iT * converter.getDeltaT() );

    // preparation for flux computations
    std::list<int> materials = { 1, 3, 4, 5 };
    IndicatorCircle3D<T> inlet( inletCenter + 2. * converter.getLatticeL() * inletNormal,
                                inletNormal, inletRadius + 2. * converter.getLatticeL() );
    IndicatorCircle3D<T> outlet0( outletCenter0 + 2. * converter.getLatticeL() * outletNormal0,
                                  outletNormal0, outletRadius0 + 2. * converter.getLatticeL() );
    IndicatorCircle3D<T> outlet1( outletCenter1 + 2. * converter.getLatticeL() * outletNormal1,
                                  outletNormal1, outletRadius1 + 2. * converter.getLatticeL() );

    // Flux of the fluid at the inlet and outlet regions
    SuperLatticePhysVelocityFlux3D<T, NSDESCRIPTOR> vFluxInflow( sLatticeNS,
        converter, superGeometry, inlet, materials );
    vFluxInflow.print( "inflow", "ml/s" );
    SuperLatticePhysPressureFlux3D<T, NSDESCRIPTOR> pFluxInflow( sLatticeNS,
        converter, superGeometry, inlet, materials );
    pFluxInflow.print( "inflow", "N", "Pa" );

    SuperLatticePhysVelocityFlux3D<T, NSDESCRIPTOR> vFluxOutflow0( sLatticeNS,
        converter, superGeometry, outlet0, materials );
    vFluxOutflow0.print( "outflow0", "ml/s" );
    SuperLatticePhysPressureFlux3D<T, NSDESCRIPTOR> pFluxOutflow0( sLatticeNS,
        converter, superGeometry, outlet0, materials );
    pFluxOutflow0.print( "outflow0", "N", "Pa" );

    SuperLatticePhysVelocityFlux3D<T, NSDESCRIPTOR> vFluxOutflow1( sLatticeNS,
        converter, superGeometry, outlet1, materials );
    vFluxOutflow1.print( "outflow1", "ml/s" );
    SuperLatticePhysPressureFlux3D<T, NSDESCRIPTOR> pFluxOutflow1( sLatticeNS,
        converter, superGeometry, outlet1, materials );
    pFluxOutflow1.print( "outflow1", "N", "Pa" );

    // Flux of particles at the inlet and outlet regions
    SuperLatticeFlux3D<T, ADDESCRIPTOR> FluxInflow( particles, superGeometry, inlet, materials );
    FluxInflow.print( "inflow" );
    SuperLatticeFlux3D<T, ADDESCRIPTOR> FluxOutflow4( particles, superGeometry, outlet0, materials );
    FluxOutflow4.print( "outflow0" );
    SuperLatticeFlux3D<T, ADDESCRIPTOR> FluxOutflow5( particles, superGeometry, outlet1, materials );
    FluxOutflow5.print( "outflow1" );

    int input[4] = {0};
    T flux3[5], flux4[5], flux5[5] = {0.};
    FluxOutflow4( flux4,input );
    FluxOutflow5( flux5,input );
    FluxInflow( flux3,input );
    // Since more diffusion is added to ensure stability the computed escaperate falls short of the real value,
    // therefore it is scaled by the factor 1.2986 computed by a simulation without drag force. This value is computed
    // for this specific setup. For further information see R.Trunk, T.Henn, W.Dörfler, H.Nirschl, M.J.Krause,
    // "Inertial Dilute Particulate Fluid Flow Simulations with an Euler-Euler Lattice Boltzmann Method"
    T escr = ( flux4[0]+flux5[0] )/flux3[0]*1.3628;
    clout << "escapeRate=" << escr << "; captureRate="<< 1-escr << std::endl;
  }
}

int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout, "main" );

  LBconverter<T> converter(
    ( int ) 3,                         // dim
    ( T ) 0.00104066/N,                // latticeL
    ( T ) 0.05/M,                      // latticeU
    ( T ) 1.5e-5,                      // charNu
    ( T ) inletRadius*2.,              // charL
    ( T ) Re*1.5e-5/( inletRadius*2 ), // charU
    ( T ) 1.225                        // charRho
  );
  converter.print();
  writeLogFile( converter, "bifurcation3d" );

  // compute relaxation parameter to solve the advection-diffusion equation in the lattice Boltzmann context
  T omegaAD = 1. / ( 4. * ( diffusion * ( converter.getDeltaT()/( converter.getDeltaX()*converter.getDeltaX() ) )
                            * ( 1./( converter.getCharU()*converter.getCharL() ) ) ) + 0.5 );

  // === 2nd Step: Prepare Geometry ===

  STLreader<T> stlReader( "../bifurcation3d.stl",converter.getLatticeL() );
  IndicatorLayer3D<T> extendedDomain( stlReader,converter.getLatticeL() );

  // Instantiation of an empty cuboidGeometry
  int noOfCuboids = std::max( 16, 4 * singleton::mpi().getSize() );
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getLatticeL(),noOfCuboids );

  // Instantiation of an empty loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  // Default instantiation of superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, NSDESCRIPTOR> sLatticeNS( superGeometry );
  SuperLattice3D<T, ADDESCRIPTOR> sLatticeAD( superGeometry );
  SuperExternal3D<T, ADDESCRIPTOR> sExternal( superGeometry, sLatticeAD, ADDESCRIPTOR<T>::ExternalField::velocityBeginsAt, ADDESCRIPTOR<T>::ExternalField::numScalars, sLatticeAD.getOverlap() );

  BGKdynamics<T, NSDESCRIPTOR> bulkDynamics( converter.getOmega(), instances::getBulkMomenta<T, NSDESCRIPTOR>() );
  ParticleAdvectionDiffusionBGKdynamics<T, ADDESCRIPTOR> bulkDynamicsAD ( omegaAD, instances::getBulkMomenta<T,ADDESCRIPTOR>() );

  sOnLatticeBoundaryCondition3D<T, NSDESCRIPTOR> sBoundaryCondition( sLatticeNS );
  createInterpBoundaryCondition3D<T, NSDESCRIPTOR>( sBoundaryCondition );

  sOnLatticeBoundaryCondition3D<T, ADDESCRIPTOR> sBoundaryConditionAD( sLatticeAD );
  createAdvectionDiffusionBoundaryCondition3D<T,ADDESCRIPTOR>( sBoundaryConditionAD );

  prepareLattice( sLatticeNS, sLatticeAD, converter, bulkDynamics, bulkDynamicsAD, sBoundaryCondition, sBoundaryConditionAD, superGeometry, omegaAD );

  AdvectionDiffusionParticleCouplingGenerator3D<T,NSDESCRIPTOR> coupling( ADDESCRIPTOR<T>::ExternalField::velocityBeginsAt );

  advDiffDragForce3D<T, NSDESCRIPTOR> dragForce( converter,radius,partRho );
  coupling.addForce( dragForce );
  sLatticeNS.addLatticeCoupling( superGeometry, 1, coupling, sLatticeAD );

  // === 4th Step: Main Loop with Timer ===
  Timer<double> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( int iT = 0; iT <= converter.numTimeSteps( maxPhysT ); ++iT ) {
    getResults( sLatticeNS, sLatticeAD, converter, iT, superGeometry, timer );
    setBoundaryValues( sLatticeNS, converter, iT, maxPhysT, superGeometry );
    sLatticeNS.executeCoupling();
    sExternal.communicate();
    sLatticeNS.collideAndStream();
    sLatticeAD.collideAndStream();
  }
  timer.stop();
  timer.printSummary();

}
