/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Patrick Nathen, Mathias J. Krause
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

#ifndef TURBULENT_F_3D_HH
#define TURBULENT_F_3D_HH

#include<vector>
#include<cmath>

#include "functors/turbulentF3D.h"
#include "functors/blockLatticeLocalF3D.h"
#include "core/superLattice3D.h"
#include "core/finiteDifference.h"
#include "core/units.h"
#include "geometry/superGeometry3D.h"
#include "utilities/vectorHelpers.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity


namespace olb {


///////////////////////////// SuperLatticeYplus3D //////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeYplus3D<T,DESCRIPTOR>::SuperLatticeYplus3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
    const LBconverter<T>& converter, SuperGeometry3D<T>& superGeometry,
    IndicatorF3D<T>& indicator, const int material )
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1),
    _superGeometry(superGeometry), _indicator(indicator), _material(material)
{
  this->getName() = "yPlus";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticeYplus3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  output[0]=T();
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<T> normalTemp(3,T());
    std::vector<T> normal(3,T());       // normalized
    T counter = T();
    T distance = T();
    if (_superGeometry.get(input) == 1) {
      for (int iPop = 1; iPop < DESCRIPTOR<T>::q; ++iPop) {
        if (_superGeometry.get(input[0],
                               input[1] + DESCRIPTOR<T>::c[iPop][0],
                               input[2] + DESCRIPTOR<T>::c[iPop][1],
                               input[3] + DESCRIPTOR<T>::c[iPop][2]) == _material) {
          counter++;
          normalTemp[0] += DESCRIPTOR<T>::c[iPop][0];
          normalTemp[1] += DESCRIPTOR<T>::c[iPop][1];
          normalTemp[2] += DESCRIPTOR<T>::c[iPop][2];
        }
      }
      if ( !util::nearZero(counter) ) {
        // get physical Coordinates at intersection

        std::vector<T> physR(3, T());
        _superGeometry.getCuboidGeometry().getPhysR(&(physR[0]), &(input[0]));

        T voxelSize = _superGeometry.getCuboidGeometry().get(globIC).getDeltaR();

        normal = util::normalize(normalTemp);

        std::vector<T> direction(normal);
        direction[0] = voxelSize*normal[0]*2.;
        direction[1] = voxelSize*normal[1]*2.;
        direction[2] = voxelSize*normal[2]*2.;

        // calculate distance to STL file
        if ( _indicator.distance(distance, physR, direction) ) {
          // call stress at this point
          T rho;
          T u[3];
          T pi[6];
          this->_sLattice.getBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix, lociy, lociz).computeRhoU(rho, u);
          this->_sLattice.getBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix, lociy, lociz).computeStress(pi);

          // Compute phys stress tau = mu*du/dx
          T omega = this->_converter.getOmega();
          T dt = this->_converter.physTime();
          T physFactor = -omega*DESCRIPTOR<T>::invCs2/rho/2./dt*this->_converter.physRho(rho)*this->_converter.getCharNu();

          //  Totel Stress projected from cell in normal direction on obstacle
          T Rx = pi[0]*physFactor*normal[0] + pi[1]*physFactor*normal[1] + pi[2]*physFactor*normal[2];
          T Ry = pi[1]*physFactor*normal[0] + pi[3]*physFactor*normal[1] + pi[4]*physFactor*normal[2];
          T Rz = pi[2]*physFactor*normal[0] + pi[4]*physFactor*normal[1] + pi[5]*physFactor*normal[2];

          // Stress appearing as pressure in corresponding direction is calculated and substracted
          T R_res_pressure = normal[0]*pi[0]*physFactor*normal[0] + normal[0]*pi[1]*physFactor*normal[1] + normal[0]*pi[2]*physFactor*normal[2]
                             +normal[1]*pi[1]*physFactor*normal[0] + normal[1]*pi[3]*physFactor*normal[1] + normal[1]*pi[4]*physFactor*normal[2]
                             +normal[2]*pi[2]*physFactor*normal[0] + normal[2]*pi[4]*physFactor*normal[1] + normal[2]*pi[5]*physFactor*normal[2];

          Rx -= R_res_pressure *normal[0];
          Ry -= R_res_pressure *normal[1];
          Rz -= R_res_pressure *normal[2];

          T tau_wall = sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
          T u_tau = sqrt(tau_wall/this->_converter.physRho(rho));
          //y_plus
          output[0] = distance*u_tau / this->_converter.getCharNu();
        } // if 4
      }
    }
  }
  return true;
}

/*template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeADM3D<T,DESCRIPTOR>::BlockLatticeADM3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, T sigma, int order, bool adaptive, const LBconverter<T>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,5), _sigma(sigma), _order(order), _adaptive(adaptive), _converter(converter)
{
  this->getName() = "ADMfilter";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  // Declaration of all variables needed for filtering

  int globX = input[0];
  int globY = input[1];
  int globZ = input[2];

  T output000[4] = {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ).computeRhoU(output000[0],output000+1);

  T outputp00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX+1, globY, globZ).computeRhoU(outputp00[0],outputp00+1);

  T output2p00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX+2, globY, globZ).computeRhoU(output2p00[0],output2p00+1);

  T outputn00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX-1, globY, globZ).computeRhoU(outputn00[0],outputn00+1);

  T output2n00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX-2, globY, globZ).computeRhoU(output2n00[0],output2n00+1);


  T output0p0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY+1, globZ).computeRhoU(output0p0[0],output0p0+1);

  T output02p0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY+2, globZ).computeRhoU(output02p0[0],output02p0+1);

  T output0n0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY-1, globZ).computeRhoU(output0n0[0],output0n0+1);

  T output02n0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY-2, globZ).computeRhoU(output02n0[0],output02n0+1);


  T output00p[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ+1).computeRhoU(output00p[0],output00p+1);

  T output002p[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ+2).computeRhoU(output002p[0],output002p+1);

  T output00n[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ-1).computeRhoU(output00n[0],output00n+1);

  T output002n[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ-2).computeRhoU(output002n[0],output002n+1);

 T relax=_sigma;

   if (_adaptive==true ){
    ///////////////////////////////////////////////DISS VERSION///////////////////////////////////////////////////

  //  T diss = dissipation(u_000)[0];

    // std::cout << "diss: "<< diss << std::endl;

  //  T* avDiss = this-> _blockLattice.get(globX, globY , globZ).getExternal(localAvDissBeginsAt);

    // // std::cout <<"avDiss:" << *avDiss << std::endl;

  //   *avDiss = (*avDiss * this->_blockLattice.getStatistics().getTime() + diss) / (this->_blockLattice.getStatistics().getTime() + (int) 1 );

    // // std::cout <<"avDiss nach mittelung:" << *avDiss << std::endl;

   //  T TKE = 0.5*(velocity(u_000)[0]*velocity(u_000)[0]+velocity(u_000)[1]*velocity(u_000)[1]+velocity(u_000)[2]+velocity(u_000)[2]);

   //  T* avTKE = this-> _blockLattice.get(globX, globY , globZ).getExternal(localAvTKEBeginsAt);

   //  *avTKE = (*avTKE * this->_blockLattice.getStatistics().getTime() + TKE) / (this->_blockLattice.getStatistics().getTime() + (int) 1 );

    // std::cout << "TKE: "<< *avTKE << std::endl;


   //  relax = sqrt((diss - *avDiss)*(diss - *avDiss));// / (*avTKE * converter.getDeltaT());
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////// Stress Version ////////////////////////////////////////////////
    T stress[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ).computeStress(stress);

    T ux = stress(u_000)[0];
    T uy = stress(u_000)[1];
    T uz = stress(u_000)[2];

    T vx = stress(u_000)[3];
    T vy = stress(u_000)[4];
    T vz = stress(u_000)[5];

    T wx = stress(u_000)[6];
    T wy = stress(u_000)[7];
    T wz = stress(u_000)[8];


    T ux = stress[0];
    T uy = stress[1];
    T uz = stress[2];

    T vx = stress[3];
    T vy = stress[4];
    T vz = stress[5];

    T wx = stress[6];
    T wy = stress[7];
    T wz = stress[8];

    T norm = sqrt(  (wx*uz + vx*uy + ux*ux) + (wy*vz + vy*vy + uy*vx) + (wz*wz + vz*wy + uz*wx) ) ;


    T* avNorm = this-> _blockLattice.get(globX, globY , globZ).getExternal(_localAvDissBeginsAt);

    *avNorm = (*avNorm * this->_blockLattice.getStatistics().getTime() + norm) / (this->_blockLattice.getStatistics().getTime() + (int) 1 );

    // relax = sigma;
    // / (*avTKE * converter.getDeltaT());


    // if (this->_blockLattice.getStatistics().getTime() >= 30000){


     relax = sqrt((norm - *avNorm)*(norm - *avNorm)) ;

      // std:: cout << "adaptive relaxation: " << relax <<  " time: "<<this->_blockLattice.getStatistics().getTime()<< endl;
    // }

   }

  if (_order==2) { // second order
    T d_0 = 6./16.;
    T d_1 = -4./16.;
    T d_2 = 1./16.;

    output[0] = output000[0] - relax*(d_2*(output2n00[0]+output02n0[0]+output002n[0]) +
                                       d_1*(outputn00[0]+output0n0[0]+output00n[0])+
                                       d_0*(output000[0]+output000[0]+output000[0])+
                                       d_1*(outputp00[0]+output0p0[0]+output00p[0])+
                                       d_2*(output2p00[0]+output02p0[0]+output002p[0]) );
    for (int i = 1; i < 4; ++i ) {
      output[i] = output000[i] - relax*(d_2*(output2n00[i] + output02n0[i] + output002n[i]) +
                                         d_1*(outputn00[i] + output0n0[i] + output00n[i])+
                                         d_0*(output000[i] + output000[i] + output000[i])+
                                         d_1*(outputp00[i] + output0p0[i] + output00p[i])+
                                         d_2*(output2p00[i] + output02p0[i] + output002p[i]) );
    }
  } else { // third order

    T output3p00[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX+3, globY, globZ).computeRhoU(output3p00[0],output3p00+1);

    T output3n00[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX-3, globY, globZ).computeRhoU(output3n00[0],output3n00+1);

    T output03p0[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY+3, globZ).computeRhoU(output03p0[0],output03p0+1);

    T output03n0[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY-3, globZ).computeRhoU(output03n0[0],output03n0+1);

    T output003p[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+3).computeRhoU(output003p[0],output003p+1);

    T output003n[4] = {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-3).computeRhoU(output003n[0],output003n+1);

    T d_0 = 5./16.;
    T d_1 = -15./64.;
    T d_2 = 3./32.;
    T d_3 = -1./64.;
    output[0] = output000[0] - _sigma*(d_3*(output3n00[0]+output03n0[0]+output003n[0])+
                                       d_2*(output2n00[0]+output02n0[0]+output002n[0]) +
                                       d_1*(outputn00[0]+output0n0[0]+output00n[0])+
                                       d_0*(output000[0]+output000[0]+output000[0])+
                                       d_1*(outputp00[0]+output0p0[0]+output00p[0])+
                                       d_2*(output2p00[0]+output02p0[0]+output002p[0])+
                                       d_3*(output3p00[0]+output03p0[0]+output003p[0]) );
    for (int i = 1; i < 4; ++i ) {
      output[i] = output000[i] - _sigma*(d_3*(output3n00[i]+output03n0[i]+output003n[i])+
                                         d_2*(output2n00[i]+output02n0[i]+output002n[i]) +
                                         d_1*(outputn00[i]+output0n0[i]+output00n[i])+
                                         d_0*(output000[i]+output000[i]+output000[i])+
                                         d_1*(outputp00[i]+output0p0[i]+output00p[i])+
                                         d_2*(output2p00[i]+output02p0[i]+output002p[i])+
                                         d_3*(output3p00[i]+output03p0[i]+output003p[i]) );
    }
  }
  output[4]=relax;
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
void BlockLatticeADM3D<T,DESCRIPTOR>::execute(const int input[])
{
  T output[5] = {T(),T(),T(),T(),T()};
  this->operator()(output,input);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineExternalField( DESCRIPTOR<T>::ExternalField::filRhoIsAt, 1, &output[0] );
  this-> _blockLattice.get(input[0],input[1],input[2]).defineExternalField( DESCRIPTOR<T>::ExternalField::localFilVelXBeginsAt, 1, &output[1]);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineExternalField( DESCRIPTOR<T>::ExternalField::localFilVelYBeginsAt, 1, &output[2]);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineExternalField( DESCRIPTOR<T>::ExternalField::localFilVelZBeginsAt, 1, &output[3]);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineExternalField( DESCRIPTOR<T>::ExternalField::localSigmaADMBeginsAt, 1, &output[4]);
}

template <typename T, template <typename U> class DESCRIPTOR>
void BlockLatticeADM3D<T,DESCRIPTOR>::execute()
{
  int nX = this-> _blockLattice.getNx();
  int nY = this-> _blockLattice.getNy();
  int nZ = this-> _blockLattice.getNz();
  int i[3];
  for (i[0]=0; i[0]<nX; ++i[0]) {
    for (i[1]=0; i[1]<nY; ++i[1]) {
      for (i[2]=0; i[2]<nZ; ++i[2]) {
        execute(i);
      }
    }
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeADM3D<T,DESCRIPTOR>::SuperLatticeADM3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, T sigma, int order, bool adaptive, const LBconverter<T>& converter)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,5), _sigma(sigma), _order(order), _adaptive(adaptive), _converter(converter)
{
  this->getName() = "ADMfilter";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticeADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    int inputLocal[3]= {};
    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), _sigma, _order,_adaptive, _converter );

    blockLatticeF(output,inputLocal);
    return true;
  } else {
    return false;
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
void SuperLatticeADM3D<T,DESCRIPTOR>::execute(SuperGeometry3D<T>& superGeometry, const int material)
{
  this->_sLattice.communicate();
  CuboidGeometry3D<T>& cGeometry =  this->_sLattice.getCuboidGeometry();
  LoadBalancer<T>& load = this->_sLattice.getLoadBalancer();
  int overlap = this->_sLattice.getOverlap();

  for (int iC = 0; iC < load.size(); ++iC) {
    BlockLatticeADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(iC), _sigma, _order, _adaptive, _converter );
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();

    int i[3];
    for (i[0]=overlap; i[0]<nX+overlap; ++i[0]) {
      for (i[1]=overlap; i[1]<nY+overlap; ++i[1]) {
        for (i[2]=overlap; i[2]<nZ+overlap; ++i[2]) {
          if (superGeometry.getExtendedBlockGeometry(iC).get(i[0],i[1],i[2]) == material) {
            blockLatticeF.execute(i);
          }
        }
      }
    }
  }
}

*/
////////////////////////BlockLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysDissipationFD3D<T,DESCRIPTOR>::BlockLatticePhysDissipationFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice ,const LBconverter<T>& _converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,3), converter(_converter)
{
  this->getName() = "PhysDissipationFD";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  //dissipation value
  T dissipation = T();

  //domain borders
  int nx = this-> _blockLattice.getNx()-1;
  int ny = this-> _blockLattice.getNy()-1;
  int nz = this-> _blockLattice.getNz()-1;

  //derivation tensor
  T dxux = T(), dxuy = T(), dxuz = T();
  T dyux = T(), dyuy = T(), dyuz = T();
  T dzux = T(), dzuy = T(), dzuz = T();

  int globX = input[0];
  int globY = input[1];
  int globZ = input[2];


  T u_000[3] = {0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ).computeU(u_000);

  //derivations in x-direction
  //boundary treatment with Second-order asymmetric gradient
  if (input[0] < 4) {
    T u_p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+1, globY, globZ).computeU(u_p00);

    T u_2p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+2, globY, globZ).computeU(u_2p00);

    dxux = fd::boundaryGradient(u_000[0], u_p00[0], u_2p00[0]);
    dxuy = fd::boundaryGradient(u_000[1], u_p00[1], u_2p00[1]);
    dxuz = fd::boundaryGradient(u_000[2], u_p00[2], u_2p00[2]);
  } else if (input[0] > nx-4) {
    T u_n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-1, globY, globZ).computeU(u_n00);

    T u_2n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-2, globY, globZ).computeU(u_2n00);

    dxux = fd::boundaryGradient(u_000[0], u_n00[0], u_2n00[0]);
    dxuy = fd::boundaryGradient(u_000[1], u_n00[1], u_2n00[1]);
    dxuz = fd::boundaryGradient(u_000[2], u_n00[2], u_2n00[2]);
  }
  //inner domain 5th order finite difference
  else {
    T u_p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+1, globY, globZ).computeU(u_p00);

    T u_2p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+2, globY, globZ).computeU(u_2p00);

    T u_3p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+3, globY, globZ).computeU(u_3p00);

    T u_4p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+4, globY, globZ).computeU(u_4p00);

    T u_n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-1, globY, globZ).computeU(u_n00);

    T u_2n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-2, globY, globZ).computeU(u_2n00);

    T u_3n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-3, globY, globZ).computeU(u_3n00);

    T u_4n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-4, globY, globZ).computeU(u_4n00);
    //TODO: Generic implementation of finite difference scheme
    dxux = ((T)672*(u_p00[0]-u_n00[0])+(T)168*(u_2n00[0]-u_2p00[0])
            +(T)32*(u_3p00[0]-u_3n00[0])+(T)3*(u_4n00[0]-u_4p00[0])) / 840;
    dxuy = ((T)672*(u_p00[1]-u_n00[1])+(T)168*(u_2n00[1]-u_2p00[1])
            +(T)32*(u_3p00[1]-u_3n00[1])+(T)3*(u_4n00[1]-u_4p00[1])) / 840;
    dxuz = ((T)672*(u_p00[2]-u_n00[2])+(T)168*(u_2n00[2]-u_2p00[2])
            +(T)32*(u_3p00[2]-u_3n00[2])+(T)3*(u_4n00[2]-u_4p00[2])) / 840;
  }
  //derivations in y-direction
  // boundary treatment with Second-order asymmetric gradient
  if (input[1] < 4) {
    T u_0p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+1, globZ).computeU(u_0p0);

    T u_02p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+2, globZ).computeU(u_02p0);

    dyux = fd::boundaryGradient(u_000[0],u_0p0[0],u_02p0[0]);
    dyuy = fd::boundaryGradient(u_000[1],u_0p0[1],u_02p0[1]);
    dyuz = fd::boundaryGradient(u_000[2],u_0p0[2],u_02p0[2]);
  } else if (input[1] > ny-4) {
    T u_0n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-1, globZ).computeU(u_0n0);

    T u_02n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-2, globZ).computeU(u_02n0);

    dyux = fd::boundaryGradient(u_000[0],u_0n0[0],u_02n0[0]);
    dyuy = fd::boundaryGradient(u_000[1],u_0n0[1],u_02n0[1]);
    dyuz = fd::boundaryGradient(u_000[2],u_0n0[2],u_02n0[2]);
  }
  //inner domain 5th order finite difference
  else {
    T u_0p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+1, globZ).computeU(u_0p0);

    T u_02p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+2, globZ).computeU(u_02p0);

    T u_03p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+3, globZ).computeU(u_03p0);

    T u_04p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+4, globZ).computeU(u_04p0);

    T u_0n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-1, globZ).computeU(u_0n0);

    T u_02n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-2, globZ).computeU(u_02n0);

    T u_03n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-3, globZ).computeU(u_03n0);

    T u_04n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-4, globZ).computeU(u_04n0);

    dyux = ((T)672*(u_0p0[0]-u_0n0[0])+(T)168*(u_02n0[0]-u_02p0[0])
            +(T)32*(u_03p0[0]-u_03n0[0])+(T)3*(u_04n0[0]-u_04p0[0])) / 840;
    dyuy = ((T)672*(u_0p0[1]-u_0n0[1])+(T)168*(u_02n0[1]-u_02p0[1])
            +(T)32*(u_03p0[1]-u_03n0[1])+(T)3*(u_04n0[1]-u_04p0[1])) / 840;
    dyuz = ((T)672*(u_0p0[2]-u_0n0[2])+(T)168*(u_02n0[2]-u_02p0[2])
            +(T)32*(u_03p0[2]-u_03n0[2])+(T)3*(u_04n0[2]-u_04p0[2])) / 840;
  }
  //derivations in z-direction
  // boundary treatment with Second-order asymmetric gradient
  if (input[2] < 4) {
    T u_00p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+1).computeU(u_00p);

    T u_002p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+2).computeU(u_002p);

    dzux = fd::boundaryGradient(u_000[0],u_00p[0],u_002p[0]);
    dzuy = fd::boundaryGradient(u_000[1],u_00p[1],u_002p[1]);
    dzuz = fd::boundaryGradient(u_000[2],u_00p[2],u_002p[2]);
  } else if (input[2] > nz-4) {
    T u_00n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-1).computeU(u_00n);

    T u_002n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-2).computeU(u_002n);

    dzux = fd::boundaryGradient(u_000[0],u_00n[0],u_002n[0]);
    dzuy = fd::boundaryGradient(u_000[1],u_00n[1],u_002n[1]);
    dzuz = fd::boundaryGradient(u_000[2],u_00n[2],u_002n[2]);
  }
  //inner domain 5th order finite difference
  else {
    T u_00p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+1).computeU(u_00p);

    T u_002p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+2).computeU(u_002p);

    T u_003p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+3).computeU(u_003p);

    T u_004p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+4).computeU(u_004p);

    T u_00n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-1).computeU(u_00n);

    T u_002n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-2).computeU(u_002n);

    T u_003n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-3).computeU(u_003n);

    T u_004n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-4).computeU(u_004n);

    dzux = ((T)672*(u_00p[0]-u_00n[0])+(T)168*(u_002n[0]-u_002p[0])
            +(T)32*(u_003p[0]-u_003n[0])+(T)3*(u_004n[0]-u_004p[0])) / 840;
    dzuy = ((T)672*(u_00p[1]-u_00n[1])+(T)168*(u_002n[1]-u_002p[1])
            +(T)32*(u_003p[1]-u_003n[1])+(T)3*(u_004n[1]-u_004p[1])) / 840;
    dzuz = ((T)672*(u_00p[2]-u_00n[2])+(T)168*(u_002n[2]-u_002p[2])
            +(T)32*(u_003p[2]-u_003n[2])+(T)3*(u_004n[2]-u_004p[2])) / 840;
  }

  //Dissipation
  dissipation =  dxux * dxux + dxuy * dxuy + dxuz * dxuz
                 +dyux * dyux + dyuy * dyuy + dyuz * dyuz
                 +dzux * dzux + dzuy * dzuy + dzuz * dzuz;

  T dt = converter.physTime();
  output[0] = dissipation * converter.getCharNu() / dt / dt;
  return true;


}

////////////////////////SuperLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDissipationFD3D<T,DESCRIPTOR>::SuperLatticePhysDissipationFD3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& _converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3),
  converter(_converter)
{
  this->getName() = "PhysDissipationFD";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    int inputLocal[3]= {};
    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticePhysDissipationFD3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), converter );

    blockLatticeF(output,inputLocal);
    return true;
  } else {
    return false;
  }


}

////////////////////////BlockLatticePhysEffectiveDissipationFD3D//////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::BlockLatticePhysEffectiveDissipationFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice ,const LBconverter<T>& _converter, T _smagoConst)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,3), converter(_converter), smagoConst(_smagoConst)
{
  this->getName() = "PhysEffectiveDissipationFD";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  //dissipation value
  T dissipation = T();

  //domain borders
  int nx = this-> _blockLattice.getNx()-1;
  int ny = this-> _blockLattice.getNy()-1;
  int nz = this-> _blockLattice.getNz()-1;

  //derivation tensor
  T dxux = T(), dxuy = T(), dxuz = T();
  T dyux = T(), dyuy = T(), dyuz = T();
  T dzux = T(), dzuy = T(), dzuz = T();

  int globX = input[0];
  int globY = input[1];
  int globZ = input[2];


  T u_000[3] = {0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ).computeU(u_000);

  //derivations in x-direction
  //boundary treatment with Second-order asymmetric gradient
  if (input[0] < 4) {
    T u_p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+1, globY, globZ).computeU(u_p00);

    T u_2p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+2, globY, globZ).computeU(u_2p00);

    dxux = fd::boundaryGradient(u_000[0], u_p00[0], u_2p00[0]);
    dxuy = fd::boundaryGradient(u_000[1], u_p00[1], u_2p00[1]);
    dxuz = fd::boundaryGradient(u_000[2], u_p00[2], u_2p00[2]);
  } else if (input[0] > nx-4) {
    T u_n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-1, globY, globZ).computeU(u_n00);

    T u_2n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-2, globY, globZ).computeU(u_2n00);

    dxux = fd::boundaryGradient(u_000[0], u_n00[0], u_2n00[0]);
    dxuy = fd::boundaryGradient(u_000[1], u_n00[1], u_2n00[1]);
    dxuz = fd::boundaryGradient(u_000[2], u_n00[2], u_2n00[2]);
  }
  //inner domain 5th order finite difference
  else {
    T u_p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+1, globY, globZ).computeU(u_p00);

    T u_2p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+2, globY, globZ).computeU(u_2p00);

    T u_3p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+3, globY, globZ).computeU(u_3p00);

    T u_4p00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX+4, globY, globZ).computeU(u_4p00);

    T u_n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-1, globY, globZ).computeU(u_n00);

    T u_2n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-2, globY, globZ).computeU(u_2n00);

    T u_3n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-3, globY, globZ).computeU(u_3n00);

    T u_4n00[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX-4, globY, globZ).computeU(u_4n00);
    //TODO: Generic implementation of finite difference scheme
    dxux = ((T)672*(u_p00[0]-u_n00[0])+(T)168*(u_2n00[0]-u_2p00[0])
            +(T)32*(u_3p00[0]-u_3n00[0])+(T)3*(u_4n00[0]-u_4p00[0])) / 840;
    dxuy = ((T)672*(u_p00[1]-u_n00[1])+(T)168*(u_2n00[1]-u_2p00[1])
            +(T)32*(u_3p00[1]-u_3n00[1])+(T)3*(u_4n00[1]-u_4p00[1])) / 840;
    dxuz = ((T)672*(u_p00[2]-u_n00[2])+(T)168*(u_2n00[2]-u_2p00[2])
            +(T)32*(u_3p00[2]-u_3n00[2])+(T)3*(u_4n00[2]-u_4p00[2])) / 840;
  }
  //derivations in y-direction
  // boundary treatment with Second-order asymmetric gradient
  if (input[1] < 4) {
    T u_0p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+1, globZ).computeU(u_0p0);

    T u_02p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+2, globZ).computeU(u_02p0);

    dyux = fd::boundaryGradient(u_000[0],u_0p0[0],u_02p0[0]);
    dyuy = fd::boundaryGradient(u_000[1],u_0p0[1],u_02p0[1]);
    dyuz = fd::boundaryGradient(u_000[2],u_0p0[2],u_02p0[2]);
  } else if (input[1] > ny-4) {
    T u_0n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-1, globZ).computeU(u_0n0);

    T u_02n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-2, globZ).computeU(u_02n0);

    dyux = fd::boundaryGradient(u_000[0],u_0n0[0],u_02n0[0]);
    dyuy = fd::boundaryGradient(u_000[1],u_0n0[1],u_02n0[1]);
    dyuz = fd::boundaryGradient(u_000[2],u_0n0[2],u_02n0[2]);
  }
  //inner domain 5th order finite difference
  else {
    T u_0p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+1, globZ).computeU(u_0p0);

    T u_02p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+2, globZ).computeU(u_02p0);

    T u_03p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+3, globZ).computeU(u_03p0);

    T u_04p0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY+4, globZ).computeU(u_04p0);

    T u_0n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-1, globZ).computeU(u_0n0);

    T u_02n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-2, globZ).computeU(u_02n0);

    T u_03n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-3, globZ).computeU(u_03n0);

    T u_04n0[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY-4, globZ).computeU(u_04n0);

    dyux = ((T)672*(u_0p0[0]-u_0n0[0])+(T)168*(u_02n0[0]-u_02p0[0])
            +(T)32*(u_03p0[0]-u_03n0[0])+(T)3*(u_04n0[0]-u_04p0[0])) / 840;
    dyuy = ((T)672*(u_0p0[1]-u_0n0[1])+(T)168*(u_02n0[1]-u_02p0[1])
            +(T)32*(u_03p0[1]-u_03n0[1])+(T)3*(u_04n0[1]-u_04p0[1])) / 840;
    dyuz = ((T)672*(u_0p0[2]-u_0n0[2])+(T)168*(u_02n0[2]-u_02p0[2])
            +(T)32*(u_03p0[2]-u_03n0[2])+(T)3*(u_04n0[2]-u_04p0[2])) / 840;
  }
  //derivations in z-direction
  // boundary treatment with Second-order asymmetric gradient
  if (input[2] < 4) {
    T u_00p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+1).computeU(u_00p);

    T u_002p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+2).computeU(u_002p);

    dzux = fd::boundaryGradient(u_000[0],u_00p[0],u_002p[0]);
    dzuy = fd::boundaryGradient(u_000[1],u_00p[1],u_002p[1]);
    dzuz = fd::boundaryGradient(u_000[2],u_00p[2],u_002p[2]);
  } else if (input[2] > nz-4) {
    T u_00n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-1).computeU(u_00n);

    T u_002n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-2).computeU(u_002n);

    dzux = fd::boundaryGradient(u_000[0],u_00n[0],u_002n[0]);
    dzuy = fd::boundaryGradient(u_000[1],u_00n[1],u_002n[1]);
    dzuz = fd::boundaryGradient(u_000[2],u_00n[2],u_002n[2]);
  }
  //inner domain 5th order finite difference
  else {
    T u_00p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+1).computeU(u_00p);

    T u_002p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+2).computeU(u_002p);

    T u_003p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+3).computeU(u_003p);

    T u_004p[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+4).computeU(u_004p);

    T u_00n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-1).computeU(u_00n);

    T u_002n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-2).computeU(u_002n);

    T u_003n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-3).computeU(u_003n);

    T u_004n[3] = {0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-4).computeU(u_004n);

    dzux = ((T)672*(u_00p[0]-u_00n[0])+(T)168*(u_002n[0]-u_002p[0])
            +(T)32*(u_003p[0]-u_003n[0])+(T)3*(u_004n[0]-u_004p[0])) / 840;
    dzuy = ((T)672*(u_00p[1]-u_00n[1])+(T)168*(u_002n[1]-u_002p[1])
            +(T)32*(u_003p[1]-u_003n[1])+(T)3*(u_004n[1]-u_004p[1])) / 840;
    dzuz = ((T)672*(u_00p[2]-u_00n[2])+(T)168*(u_002n[2]-u_002p[2])
            +(T)32*(u_003p[2]-u_003n[2])+(T)3*(u_004n[2]-u_004p[2])) / 840;
  }

  //Dissipation
  dissipation =  dxux * dxux + dxuy * dxuy + dxuz * dxuz
                 +dyux * dyux + dyuy * dyuy + dyuz * dyuz
                 +dzux * dzux + dzuy * dzuy + dzuz * dzuz;

  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. / converter.getOmega();
  /// Turbulent realaxation time
  T nuLattice = converter.getLatticeNu();
  T preFactor = (smagoConst*smagoConst)*DESCRIPTOR<T>::invCs2*DESCRIPTOR<T>::invCs2*2*sqrt(2);
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor/rho*PiNeqNorm) - tau_mol);
  T nuTurbulentLattice = tau_turb/DESCRIPTOR<T>::invCs2;

  T dt = converter.physTime();
  output[0] = dissipation * (nuTurbulentLattice + nuLattice)* converter.getCharNu() / converter.getLatticeNu() / dt / dt;

  return true;


}

////////////////////////SuperLatticePhysEffectiveDissipationFD3D//////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::SuperLatticePhysEffectiveDissipationFD3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& _converter, T _smagoConst) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3),
  converter(_converter), smagoConst(_smagoConst)
{
  this->getName() = "PhysEffectiveDissipationFD";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    int inputLocal[3]= {};
    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), converter, smagoConst );

    blockLatticeF(output,inputLocal);
    return true;
  } else {
    return false;
  }


}
/*
////////////////////////BlockLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeSigmaADM3D<T,DESCRIPTOR>::BlockLatticeSigmaADM3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice )
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,3)
{
  this->getName() = "SigmaADM3D";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticeSigmaADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

 T* sigma = this-> _blockLattice.get(input[0], input[1] , input[2]).getExternal(_localSigmaADMBeginsAt);
   output[0] = *sigma;

  return true;

}


////////////////////////SuperLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeSigmaADM3D<T,DESCRIPTOR>::SuperLatticeSigmaADM3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3)
{
  this->getName() = "SigmaADM3D";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticeSigmaADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    int inputLocal[3]= {};
    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeSigmaADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)));

    blockLatticeF(output,inputLocal);
    return true;
  } else {
    return false;
  }


}

*/

} // end namespace olb
#endif


