/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "dynamics.h"
#include "dynamics.hh"
#include "latticeDescriptors.h"
#include "latticeDescriptors.hh"

namespace olb {

// Efficient specialization for D3Q19 lattice
template<>
void BGKdynamics<double, descriptors::D3Q19Descriptor >::collide (
  Cell<double,descriptors::D3Q19Descriptor >& cell,
  LatticeStatistics<double>& statistics )
{
  typedef double T;
  typedef descriptors::D3Q19DescriptorBase<T> L;

  T rho, u[3];
  this->_momenta.computeRhoU(cell, rho, u);

  T one_m_omega = (T)1 - _omega;
  T t0_omega = L::t[0]*_omega; // weight for i=0
  T t1_omega = L::t[1]*_omega; // weight for i=1,2,3,10,11,12
  T t4_omega = L::t[4]*_omega; // weight for i=4,5,6,7,8,9,13,14,15,16,17,18

  T uSqr     = u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; // compute of usqr

  T u3x      = (T)3*u[0]; // compute of 3*ux
  T u3y      = (T)3*u[1]; // compute of 3*uy
  T u3z      = (T)3*u[2]; // compute of 3*uz

  T u3xSqr_  = .5*u3x*u3x; // compute 9/2*(ux)²
  T u3ySqr_  = .5*u3y*u3y; // compute 9/2*(uy)²
  T u3zSqr_  = .5*u3z*u3z;  // compute 9/2*(uz)²

  T u3xu3y_  = u3x*u3y; // compute 9(ux*uy)
  T u3xu3z_  = u3x*u3z; // compute 9(ux*uz)
  T u3yu3z_  = u3y*u3z; // compute 9(uy*uz)

  T C1 = (T)1 + (T)3*uSqr; // compute 1+3*((ux)²+(uy)²+(uz)²)
  T C2, C3;

  //**************************case i=0
  C3 = -u3xSqr_ - u3ySqr_ - u3zSqr_; // compute -9/2*((ux)²+(uy)²+(uz)²)

  // compute of feq0=rho*t0*(c1+c3)=rho*t0*(1-3/2*((ux)²+(uy)²+(uz)²))
  cell[0] *= one_m_omega;
  cell[0] += t0_omega*(rho*(C1 + C3) - (T)1);

  //**************************case i=1 and i=10
  C2 = -u3x; // compute -3*ux
  C3 = -u3ySqr_ - u3zSqr_; // compute -9/2*((uy)²+(uz)²)

  // compute of feq1=rho*t1*(c1+c2+c3)=rho*t0*(1-3*ux+9/2*(ux)²-3/2*((ux)²+(uy)²+(uz)²))
  cell[1]  *= one_m_omega;
  cell[1]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq10=rho*t1*(c1-c2+c3)=rho*t0*(1+3*ux+9/2*(ux)²-3/2*((ux)²+(uy)²+(uz)²))
  cell[10] *= one_m_omega;
  cell[10] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=2 and i=11
  C2 = -u3y; // compute -3*uy
  C3 = -u3xSqr_ - u3zSqr_; // compute -9/2*((ux)²+(uz)²)

  // compute of feq2=rho*t1*(c1-c2+c3)=rho*t1*(1-3*uy+9/2*(uy)²-3/2*((ux)²+(uy)²+(uz)²))
  cell[2]  *= one_m_omega;
  cell[2]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq2=rho*t1*(c1+c2+c3)=rho*t1*(1+3*uy+9/2*(uy)²-3/2*((ux)²+(uy)²+(uz)²))
  cell[11] *= one_m_omega;
  cell[11] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=3 and i=12
  C2 = -u3z; // compute -3*uz
  C3 = -u3xSqr_ - u3ySqr_; // compute -9/2*((ux)²+(uy)²)

  // compute of feq3=rho*t1*(c1+c2+c3)=rho*t1*(1-3*uz+9/2*(uz)²-3/2*((ux)²+(uy)²+(uz)²))
  cell[3]  *= one_m_omega;
  cell[3]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq12=rho*t1*(c1-c2+c3)=rho*t1*(1+3*uz+9/2*(uz)²-3/2*((ux)²+(uy)²+(uz)²))
  cell[12] *= one_m_omega;
  cell[12] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=4 and i=13
  C2 = -u3x - u3y; // compute -3*(uz+uy)
  C3 = u3xu3y_ - u3zSqr_; // compute 9*(ux*uy)-9/2*(uz)²

  // compute of feq4=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux+uy)+9/2*((ux)²+(uy)²+2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
  cell[4]  *= one_m_omega;
  cell[4]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq13=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux+uy)+9/2*((ux)²+(uy)²+2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
  cell[13] *= one_m_omega;
  cell[13] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=5 and i=14
  C2 = -u3x + u3y; // compute -3*(ux+uy)
  C3 = -u3xu3y_ - u3zSqr_; // compute -9*(ux*uy)-9/2*(uz)²

  // compute of feq5=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux-uy)+9/2*((ux)²+(uy)²-2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
  cell[5]  *= one_m_omega;
  cell[5]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq14=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux-uy)+9/2*((ux)²+(uy)²-2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
  cell[14] *= one_m_omega;
  cell[14] +=t4_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=6 and i=15
  C2 = -u3x - u3z; // compute -3*(ux+uz)
  C3 = u3xu3z_ - u3ySqr_; // compute 9*(ux*uz)-9/2*(uy)²

  // compute of feq6=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux+uz)+9/2*((ux)²+(uy)²+2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[6]  *= one_m_omega;
  cell[6]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq15=rho*t4*(c1-c2+c3)=rho*t15*(1+3*(ux+uz)+9/2*((ux)²+(uy)²+2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[15] *= one_m_omega;
  cell[15] +=t4_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=7 and i=16
  C2 = -u3x + u3z; // compute -3*(ux-uz)
  C3 = -u3xu3z_ - u3ySqr_; // compute -9*(ux*uz)-9/2*(uy)²

  // compute of feq7=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux-uz)+9/2*((ux)²+(uz)²-2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[7]  *= one_m_omega;
  cell[7]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq16=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux-uz)+9/2*((ux)²+(uz)²-2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[16] *= one_m_omega;
  cell[16] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=7 and i=16
  C2 = -u3y - u3z; // compute -3*(uy+uz)
  C3 = u3yu3z_ - u3xSqr_; // compute 9*(uy*uz)-9/2*(ux)²

  // compute of feq8=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(uy+uz)+9/2*((uz)²+(uy)²+2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[8]  *= one_m_omega;
  cell[8]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq17=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(uy+uz)+9/2*((uz)²+(uy)²+2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[17] *= one_m_omega;
  cell[17] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

  //************************case i=7 and i=16
  C2 = -u3y + u3z; // compute -3*(uy-uz)
  C3 = -u3yu3z_ - u3xSqr_; // compute -9*(uy*uz)-9/2*(ux)²

  // compute of feq9=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(uy-uz)+9/2*((uy)²+(uz)²-2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[9]  *= one_m_omega;
  cell[9]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

  // compute of feq18=rho*t4*(c1+c2+c3)=rho*t4*(1+3*(uy-uz)+9/2*((uy)²+(uz)²-2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
  cell[18] *= one_m_omega;
  cell[18] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

  if (cell.takesStatistics()) {
    statistics.incrementStats(rho,uSqr);
  }
};

template struct Dynamics<double, descriptors::D2Q9Descriptor>;
template struct Momenta<double, descriptors::D2Q9Descriptor>;
template class BasicDynamics<double, descriptors::D2Q9Descriptor>;
template class BGKdynamics<double, descriptors::D2Q9Descriptor>;
template class ConstRhoBGKdynamics<double, descriptors::D2Q9Descriptor>;
template class IncBGKdynamics<double, descriptors::D2Q9Descriptor>;
template class RLBdynamics<double, descriptors::D2Q9Descriptor>;
template class CombinedRLBdynamics<double, descriptors::D2Q9Descriptor,
                                   RLBdynamics<double, descriptors::D2Q9Descriptor> >;
template class CombinedRLBdynamics<double, descriptors::D2Q9Descriptor,
                                   BGKdynamics<double, descriptors::D2Q9Descriptor> >;
template class CombinedRLBdynamics<double, descriptors::D2Q9Descriptor,
                                   ConstRhoBGKdynamics<double, descriptors::D2Q9Descriptor> >;
template struct BulkMomenta<double, descriptors::D2Q9Descriptor>;
template class BounceBack<double, descriptors::D2Q9Descriptor>;
template class BounceBackVelocity<double, descriptors::D2Q9Descriptor>;
template class BounceBackAnti<double, descriptors::D2Q9Descriptor>;
template class NoDynamics<double, descriptors::D2Q9Descriptor>;
template class OffDynamics<double, descriptors::D2Q9Descriptor>;
template class ZeroDistributionDynamics<double, descriptors::D2Q9Descriptor>;

template struct Dynamics<double, descriptors::D3Q19Descriptor>;
template struct Momenta<double, descriptors::D3Q19Descriptor>;
template class BasicDynamics<double, descriptors::D3Q19Descriptor>;
template class BGKdynamics<double, descriptors::D3Q19Descriptor>;
template class ConstRhoBGKdynamics<double, descriptors::D3Q19Descriptor>;
template class IncBGKdynamics<double, descriptors::D3Q19Descriptor>;
template class RLBdynamics<double, descriptors::D3Q19Descriptor>;
template class CombinedRLBdynamics<double, descriptors::D3Q19Descriptor,
                                   RLBdynamics<double, descriptors::D3Q19Descriptor> >;
template class CombinedRLBdynamics<double, descriptors::D3Q19Descriptor,
                                   BGKdynamics<double, descriptors::D3Q19Descriptor> >;
template class CombinedRLBdynamics<double, descriptors::D3Q19Descriptor,
                                   ConstRhoBGKdynamics<double, descriptors::D3Q19Descriptor> >;
template struct BulkMomenta<double, descriptors::D3Q19Descriptor>;
template class BounceBack<double, descriptors::D3Q19Descriptor>;
template class BounceBackVelocity<double, descriptors::D3Q19Descriptor>;
template class BounceBackAnti<double, descriptors::D3Q19Descriptor>;
template class NoDynamics<double, descriptors::D3Q19Descriptor>;
template class OffDynamics<double, descriptors::D3Q19Descriptor>;
template class ZeroDistributionDynamics<double, descriptors::D3Q19Descriptor>;


namespace instances {

template BulkMomenta<double, descriptors::D2Q9Descriptor>& getBulkMomenta();

template BounceBack<double, descriptors::D2Q9Descriptor>& getBounceBack();

template BounceBackAnti<double, descriptors::D2Q9Descriptor>& getBounceBackAnti(const double rho);

template BounceBackVelocity<double, descriptors::D2Q9Descriptor>& getBounceBackVelocity(const double rho, const double u[2]);

template NoDynamics<double, descriptors::D2Q9Descriptor>& getNoDynamics(double rho);

template ZeroDistributionDynamics<double, descriptors::D2Q9Descriptor>& getZeroDistributionDynamics();

template BulkMomenta<double, descriptors::D3Q19Descriptor>& getBulkMomenta();

template BounceBack<double, descriptors::D3Q19Descriptor>& getBounceBack();

template BounceBackVelocity<double, descriptors::D3Q19Descriptor>& getBounceBackVelocity(const double rho, const double u[3]);

template BounceBackAnti<double, descriptors::D3Q19Descriptor>& getBounceBackAnti(const double rho);

template NoDynamics<double, descriptors::D3Q19Descriptor>& getNoDynamics(double rho);

template ZeroDistributionDynamics<double, descriptors::D3Q19Descriptor>& getZeroDistributionDynamics();
}

}
