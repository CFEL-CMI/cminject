/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
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

/* rayleighBenard2d.cpp:
 * Rayleigh-Benard convection rolls in 2D, simulated with
 * the thermal LB model by Z. Guo e.a., between a hot plate at
 * the bottom and a cold plate at the top.
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

#define NSDESCRIPTOR ForcedD2Q9Descriptor
#define TDESCRIPTOR AdvectionDiffusionD2Q5Descriptor

const int maxIter  = 1000000000;
const int saveIter = 2500;

void prepareGeometry( SuperGeometry2D<T>& superGeometry,
                      AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary

  superGeometry.rename( 0,1 );

  Vector<T,2> origin1( -2., T() );
  Vector<T,2> origin2( -2., converter.getLy()*converter.getN() );
  Vector<T,2> origin3( converter.getLx()*converter.getN()/2., 1. );
  Vector<T,2> extend1( converter.getLx()*converter.getN()+4, T() );
  Vector<T,2> extend2( 1 );

  IndicatorCuboid2D<T> bottom( extend1, origin1 );
  IndicatorCuboid2D<T> top( extend1, origin2 );
  IndicatorCuboid2D<T> perturbation( extend2, origin3 );

  superGeometry.rename( 1,2,bottom );
  superGeometry.rename( 1,3,top );
  superGeometry.rename( 1,4,perturbation );

  // Removes all not needed boundary voxels outside the surface
  //superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter,
                     SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                     ForcedBGKdynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>& NSboundaryCondition,
                     sOnLatticeBoundaryCondition2D<T,TDESCRIPTOR>& TboundaryCondition,
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );

  double Tomega  = converter.getOmegaT();
  double NSomega = converter.getOmegaNS();

  // define lattice Dynamics
  clout << "defining dynamics" << endl;

  ADlattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>() );
  NSlattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>() );

  ADlattice.defineDynamics( superGeometry, 1, &advectionDiffusionBulkDynamics );
  ADlattice.defineDynamics( superGeometry, 2, &advectionDiffusionBulkDynamics );
  ADlattice.defineDynamics( superGeometry, 3, &advectionDiffusionBulkDynamics );
  ADlattice.defineDynamics( superGeometry, 4, &advectionDiffusionBulkDynamics );
  NSlattice.defineDynamics( superGeometry, 1, &bulkDynamics );
  NSlattice.defineDynamics( superGeometry, 2, &bulkDynamics );
  NSlattice.defineDynamics( superGeometry, 3, &bulkDynamics );
  NSlattice.defineDynamics( superGeometry, 4, &bulkDynamics );

  // sets boundary
  NSboundaryCondition.addVelocityBoundary( superGeometry, 2, NSomega );
  NSboundaryCondition.addVelocityBoundary( superGeometry, 3, NSomega );
  TboundaryCondition.addTemperatureBoundary( superGeometry, 2, Tomega );
  TboundaryCondition.addTemperatureBoundary( superGeometry, 3, Tomega );
}

void setBoundaryValues( AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter,
                        SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                        SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                        int iT, SuperGeometry2D<T>& superGeometry ) {

  if ( iT == 0 ) {

    typedef advectionDiffusionLbHelpers<T,TDESCRIPTOR> TlbH;

    // for each material set the defineRhoU and the Equilibrium
    std::vector<T> zero( 2,T() );
    T zerovel[TDESCRIPTOR<T>::d] = {0.,0.};
    AnalyticalConst2D<T,T> u( zero );
    AnalyticalConst2D<T,T> rho( 1. );
    AnalyticalConst2D<T,T> force( zero );
    NSlattice.defineRhoU( superGeometry, 1, rho, u );
    NSlattice.iniEquilibrium( superGeometry, 1, rho, u );
    NSlattice.defineExternalField( superGeometry, 1,
                                   NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                   NSDESCRIPTOR<T>::ExternalField::sizeOfForce, force );
    NSlattice.defineRhoU( superGeometry, 2, rho, u );
    NSlattice.iniEquilibrium( superGeometry, 2, rho, u );
    NSlattice.defineExternalField( superGeometry, 2,
                                   NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                   NSDESCRIPTOR<T>::ExternalField::sizeOfForce, force );
    NSlattice.defineRhoU( superGeometry, 3, rho, u );
    NSlattice.iniEquilibrium( superGeometry, 3, rho, u );
    NSlattice.defineExternalField( superGeometry, 3,
                                   NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                   NSDESCRIPTOR<T>::ExternalField::sizeOfForce, force );
    NSlattice.defineRhoU( superGeometry, 4, rho, u );
    NSlattice.iniEquilibrium( superGeometry, 4, rho, u );
    NSlattice.defineExternalField( superGeometry, 4,
                                   NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                   NSDESCRIPTOR<T>::ExternalField::sizeOfForce, force );

    T Tcold    = converter.getTcold();
    T Thot     = converter.getThot();
    T Tperturb = Thot+Thot/( T )5;

    AnalyticalConst2D<T,T> Cold( Tcold );
    AnalyticalConst2D<T,T> Hot( Thot );
    AnalyticalConst2D<T,T> Perturb( Tperturb );

    std::vector<T> tEqCold( TDESCRIPTOR<T>::q );
    std::vector<T> tEqHot( TDESCRIPTOR<T>::q );
    std::vector<T> tEqPerturb( TDESCRIPTOR<T>::q );

    for ( int iPop = 0; iPop < TDESCRIPTOR<T>::q; ++iPop ) {
      tEqCold[iPop] = TlbH::equilibrium( iPop,Tcold,zerovel );
      tEqHot[iPop] = TlbH::equilibrium( iPop,Thot,zerovel );
      tEqPerturb[iPop] = TlbH::equilibrium( iPop,Tperturb,zerovel );
    }

    AnalyticalConst2D<T,T> EqCold( tEqCold );
    AnalyticalConst2D<T,T> EqHot( tEqHot );
    AnalyticalConst2D<T,T> EqPerturb( tEqPerturb );

    ADlattice.defineRho( superGeometry, 1, Cold );
    ADlattice.definePopulations( superGeometry, 1, EqCold );
    ADlattice.defineExternalField( superGeometry, 1,
                                   TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
                                   TDESCRIPTOR<T>::ExternalField::sizeOfVelocity, u );
    ADlattice.defineRho( superGeometry, 2, Hot );
    ADlattice.definePopulations( superGeometry, 2, EqHot );
    ADlattice.defineExternalField( superGeometry, 2,
                                   TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
                                   TDESCRIPTOR<T>::ExternalField::sizeOfVelocity, u );
    ADlattice.defineRho( superGeometry, 3, Cold );
    ADlattice.definePopulations( superGeometry, 3, EqCold );
    ADlattice.defineExternalField( superGeometry, 3,
                                   TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
                                   TDESCRIPTOR<T>::ExternalField::sizeOfVelocity, u );
    ADlattice.defineRho( superGeometry, 4, Perturb );
    ADlattice.definePopulations( superGeometry, 4, EqPerturb );
    ADlattice.defineExternalField( superGeometry, 4,
                                   TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
                                   TDESCRIPTOR<T>::ExternalField::sizeOfVelocity, u );

    // Make the lattice ready for simulation
    NSlattice.initialize();
    ADlattice.initialize();
  }
}

void getResults( AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter,
                 SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice2D<T, TDESCRIPTOR>& ADlattice, int iT,
                 SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriterNS( "rayleighBenard2dNSlattice" );
  SuperVTMwriter2D<T> vtmWriterAD( "rayleighBenard2dADlattice" );

  if ( iT == 0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid( NSlattice );
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank( NSlattice );
    vtmWriterNS.write( geometry );
    vtmWriterNS.write( cuboid );
    vtmWriterNS.write( rank );
    vtmWriterNS.createMasterFile();
    vtmWriterAD.createMasterFile();
  }

  // Writes the vtm files
  if ( iT%saveIter == 0 ) {
    clout << "Writing stats at time " << iT << "." << endl;
    clout << ADlattice.getStatistics().getAverageEnergy() << endl;
    clout << "Writing vtm..." << endl;
    SuperLatticeVelocity2D<T, NSDESCRIPTOR> velocity( NSlattice );
    SuperLatticeDensity2D<T, TDESCRIPTOR> density( ADlattice );
    vtmWriterNS.addFunctor( velocity );
    vtmWriterAD.addFunctor( density );
    vtmWriterNS.write( iT );
    vtmWriterAD.write( iT );

    BlockLatticeReduction2D<T, TDESCRIPTOR> planeReduction( density );
    BlockGifWriter<T> gifWriter;
    gifWriter.write( planeReduction, iT, "density" );
  }

}

T computeNusselt( SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                  SuperLattice2D<T,TDESCRIPTOR> &ADlattice,
                  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> &converter ) {
  T u_T = T();
  for ( int iC = 0; iC < ADlattice.getLoadBalancer().size(); iC++ ) {
    int nx = ADlattice.getBlockLattice( iC ).getNx();
    int ny = ADlattice.getBlockLattice( iC ).getNy();

    for ( int iX = 0; iX < nx; ++iX ) {
      for ( int iY = 0; iY < ny; ++iY ) {
        T uy[2];
        NSlattice.getBlockLattice( iC ).get( iX,iY ).computeU( uy );
        u_T += uy[1] * ADlattice.getBlockLattice( iC ).get( iX,iY ).computeRho();
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast( u_T, MPI_SUM );
#endif
  T nusselt = ( T )1 + u_T*converter.getDeltaX() / ( converter.getKappa() * ( converter.getThot()-converter.getTcold() ) );

  return nusselt;
}

int main( int argc, char *argv[] ) {

  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"main" );
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );

  double Ra    = 10000;
  double Pr    = 1;
  int N        = 50;
  double dt    = 0.001;

  if ( argc == 5 ) {
    Ra    = atof( argv[1] );
    Pr    = atof( argv[2] );
    N     = atoi( argv[3] );
    dt    = atof( argv[4] );
  } else {
    clout << "Warning : Wrong number of parameters specified." << endl;
    clout << "1 : Rayleigh.    default: " << Ra << endl;
    clout << "2 : Prandtl.     default: " << Pr << endl;
    clout << "3 : N.           default: " << N << endl;
    clout << "4 : Delta t.     default: " << dt << endl;
    clout << "Now using default values instead.";
  }

  clout << "Starting simulation with parameters Ra=" << Ra << " Pr=" << Pr << " N=" << N << " dt=" << dt << endl;

  AdvectionDiffusionUnitLB<T,NSDESCRIPTOR,TDESCRIPTOR> converter(
    Ra,   // Ra
    Pr,   // Pr
    0.0,  // Tcold
    1.0,  // Thot
    N,    // N
    dt,   // dt
    2.0,  // lx
    1.0   // ly
  );
  writeLogFile<T,NSDESCRIPTOR,TDESCRIPTOR>( converter,"rayleighBenard2d" );

  const double Raprova = converter.getN() * converter.getN() *
                         converter.getN() * converter.getDeltaTemperature() *
                         converter.getGravity() / ( converter.getNu()*converter.getKappa() );

  const double Prprova = converter.getNu() / converter.getKappa();

  clout << Raprova << " " << Prprova << endl;
  clout << converter.getOmegaNS() << " " << converter.getOmegaT() << endl;

  int nx = converter.getNx();
  int ny = converter.getNy();

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cGeometry( 0, 0, 1, nx, ny, singleton::mpi().getSize() );
#else
  CuboidGeometry2D<T> cGeometry( 0, 0, 1, nx, ny, 1 );
#endif

  cGeometry.setPeriodicity( true, false );

  HeuristicLoadBalancer<T> loadBalancer( cGeometry );

  SuperGeometry2D<T> superGeometry( cGeometry, loadBalancer, 2 );

  prepareGeometry( superGeometry, converter );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice2D<T, TDESCRIPTOR> ADlattice( superGeometry );
  SuperLattice2D<T, NSDESCRIPTOR> NSlattice( superGeometry );

  sOnLatticeBoundaryCondition2D<T, NSDESCRIPTOR> NSboundaryCondition( NSlattice );
  createLocalBoundaryCondition2D<T, NSDESCRIPTOR>( NSboundaryCondition );

  sOnLatticeBoundaryCondition2D<T, TDESCRIPTOR> TboundaryCondition( ADlattice );
  createAdvectionDiffusionBoundaryCondition2D<T, TDESCRIPTOR>( TboundaryCondition );

  ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getOmegaNS(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>() );

  AdvectionDiffusionRLBdynamics<T, TDESCRIPTOR> TbulkDynamics (
    converter.getOmegaT(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>()
  );

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<T> dir;
  dir.push_back( T() );
  dir.push_back( ( T )1 );

  NavierStokesAdvectionDiffusionCouplingGenerator2D<T,NSDESCRIPTOR>
  coupling( 0, nx-1, 0, ny-1, converter.getGravity(),
            converter.getT0(), converter.getDeltaTemperature(), dir );

  NSlattice.addLatticeCoupling( superGeometry, 1, coupling, ADlattice );
  NSlattice.addLatticeCoupling( superGeometry, 2, coupling, ADlattice );
  NSlattice.addLatticeCoupling( superGeometry, 3, coupling, ADlattice );
  NSlattice.addLatticeCoupling( superGeometry, 4, coupling, ADlattice );

  prepareLattice( converter,
                  NSlattice, ADlattice,
                  NSbulkDynamics, TbulkDynamics,
                  NSboundaryCondition, TboundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===

  util::ValueTracer<T> converge( 0.01,( T )ny,1.0e-5 );
  int iT = 0;
  clout << "starting simulation" << endl;
  for ( iT = 0; iT < maxIter; ++iT ) {

    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, NSlattice, ADlattice, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    ADlattice.collideAndStream();
    NSlattice.collideAndStream();

    converge.takeValue( ADlattice.getStatistics().getAverageEnergy(),true );

    NSlattice.executeCoupling();
    //ADlattice.executeCoupling();

    // === 7th Step: Computation and Output of the Results ===
    getResults( converter, NSlattice, ADlattice, iT, superGeometry );
  }

  clout << "Time " << iT << "." << std::endl;
  clout << "Nusselt = " << computeNusselt( NSlattice, ADlattice, converter ) << endl;
}
