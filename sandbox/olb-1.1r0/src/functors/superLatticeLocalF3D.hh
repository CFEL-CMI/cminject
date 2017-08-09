/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2014 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef SUPER_LATTICE_LOCAL_F_3D_HH
#define SUPER_LATTICE_LOCAL_F_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm

#include "functors/superBaseF3D.h"
#include "functors/superLatticeLocalF3D.h"
#include "functors/indicator/indicatorBaseF3D.hh"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry3D.h"

namespace olb {



template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFpop3D<T, DESCRIPTOR>::SuperLatticeFpop3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR<T>::q)
{
  this->getName() = "fPop";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticeFpop3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFpop3D<T, DESCRIPTOR>::operator()( T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeDissipation3D<T, DESCRIPTOR>::SuperLatticeDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "dissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticeDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeDissipation3D<T, DESCRIPTOR>::operator()( T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysDissipation3D<T, DESCRIPTOR>::SuperLatticePhysDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticePhysDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDissipation3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::SuperLatticeEffevtiveDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter, T smagoConst)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "EffevtiveDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter, smagoConst));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::operator()( T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::SuperLatticePhysEffevtiveDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter, T smagoConst)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physEffevtiveDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter, smagoConst) );
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeDensity3D<T, DESCRIPTOR>::SuperLatticeDensity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "density";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticeDensity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeDensity3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeVelocity3D<T, DESCRIPTOR>::SuperLatticeVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "velocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticeVelocity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeVelocity3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeStrainRate3D<T, DESCRIPTOR>::SuperLatticeStrainRate3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 9), _converter(converter)
{
  this->getName() = "strainRate";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticeStrainRate3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeStrainRate3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysStrainRate3D<T, DESCRIPTOR>::SuperLatticePhysStrainRate3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 9)
{
  this->getName() = "physStrainRate";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticePhysStrainRate3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysStrainRate3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeGeometry3D<T, DESCRIPTOR>::SuperLatticeGeometry3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "geometry";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new  BlockLatticeGeometry3D<T, DESCRIPTOR>(
                              this->_sLattice.getBlockLattice(iC),
                              this->_superGeometry.getBlockGeometry(iC),
                              _material) );
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeGeometry3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeRank3D<T, DESCRIPTOR>::SuperLatticeRank3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "rank";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back( new BlockLatticeRank3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)) );
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeRank3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeCuboid3D<T, DESCRIPTOR>::SuperLatticeCuboid3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "cuboid";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back( new BlockLatticeCuboid3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                             this->_sLattice.getLoadBalancer().glob(iC)) );
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeCuboid3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressure3D<T, DESCRIPTOR>::SuperLatticePhysPressure3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physPressure";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticePhysPressure3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysPressure3D<T, DESCRIPTOR>::operator()( T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocity3D<T, DESCRIPTOR>::SuperLatticePhysVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter, bool print)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3), _print(print)
{
  this->getName() = "physVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticePhysVelocity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                            this->_converter, _print));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeExternal3D<T, DESCRIPTOR>::SuperLatticeExternal3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, int start, int size)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, size)
{
  this->getName() = "ExtField";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); iC++ ) {
    this->_blockF.push_back(new BlockLatticeExternal3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), start, size));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeExternal3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysExternal3D<T, DESCRIPTOR>::SuperLatticePhysExternal3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3)
{
  this->getName() = "physExtField";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticePhysExternal3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysExternal3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysExternalPorosity3D<T,DESCRIPTOR>::SuperLatticePhysExternalPorosity3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1)
{
  this->getName() = "ExtPorosityField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  int globIC = input[0];
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    int inputLocal[3] = { };

    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticePhysExternalPorosity3D<T, DESCRIPTOR> blockLatticeF(
      this->_sLattice.getExtendedBlockLattice(
        this->_sLattice.getLoadBalancer().loc(globIC)),
      this->_converter);
    blockLatticeF(output, inputLocal);
//    std::cout << output[0] << std::endl;
//    std::cout << output[0];


    return true;
  } else {
    return false;
  }
  return false;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysExternalVelocity3D<T, DESCRIPTOR>::SuperLatticePhysExternalVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtVelocityField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysExternalVelocity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  int globIC = input[0];
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    int inputLocal[3] = { };

    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    BlockLatticePhysExternalVelocity3D<T, DESCRIPTOR> blockLatticeF(
      this->_sLattice.getExtendedBlockLattice(
        this->_sLattice.getLoadBalancer().loc(globIC)),
      this->_converter);
    blockLatticeF(output, inputLocal);
    return true;
  } else {
    return false;
  }
  return false;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>::SuperLatticePhysExternalParticleVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtPartVelField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  int globIC = input[0];
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    int inputLocal[3] = { };

    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    BlockLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR> blockLatticeF(
      this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)),
      this->_converter);
    blockLatticeF(output, inputLocal);
    return true;
  } else {
    return false;
  }
  return false;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physBoundaryForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(
      new BlockLatticePhysBoundaryForce3D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC), _superGeometry.getBlockGeometry(iC), _material, this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysBoundaryForce3D<T, DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


//TODO
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysBoundaryForceIndicator3D<T,DESCRIPTOR>::SuperLatticePhysBoundaryForceIndicator3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
 ParticleIndicatorF3D<T,T>& indicator, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physBoundaryForceIndicator";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysBoundaryForceIndicator3D<T,DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];
  T inside[1];

  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;

    std::vector<T> posIn(3,T());
    posIn = this->_superGeometry.getPhysR(globIC, locix, lociy, lociz);
    _indicator(inside, &(posIn[0]) );
    if ( !util::nearZero(inside[0]) ) {
      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
        // Get direction
        const int* c = DESCRIPTOR<T>::c[iPop];
        // Get next cell located in the current direction
        // Check if the next cell is a fluid node
//        std::vector<int> input2(input,input+4);
//        input2[1] = input[1] + c[0];
//        input2[2] = input[2] + c[1];
//        input2[3] = input[3] + c[2];
//
//        std::vector<T> posOut(3,T());
//        posOut = this->_superGeometry.getPhysR(globIC, input2[1], input2[2], input2[3]);
//        _indicator(inside, &(posOut[0]) );
//        if ( inside[0] == 0) {
        int overlap = this->_sLattice.getOverlap();
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
        // Update force
        output[0] -= c[0]*f;
        output[1] -= c[1]*f;
        output[2] -= c[2]*f;
//        }
      }
      output[0] = this->_converter.physForce(output[0]);
      output[1] = this->_converter.physForce(output[1]);
      output[2] = this->_converter.physForce(output[2]);
//      return true;
    }
    //else {
    return true;
  } else {
    return false;
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysBoundaryTorqueIndicator3D<T,DESCRIPTOR>::SuperLatticePhysBoundaryTorqueIndicator3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
 ParticleIndicatorF3D<T,T>& indicator, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physBoundaryTorqueIndicator";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysBoundaryTorqueIndicator3D<T,DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];
  T inside[1];

  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<T> force(3, T());
    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;
    std::vector<T> posIn(3,T());
    posIn = this->_superGeometry.getPhysR(globIC, locix, lociy, lociz);
    _indicator(inside, &(posIn[0]) );
    if ( !util::nearZero(inside[0]) ) {
      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
        // Get direction
        const int* c = DESCRIPTOR<T>::c[iPop];
        // Get next cell located in the current direction
        // Check if the next cell is a fluid node
        std::vector<int> input2(input,input+4);
        input2[1] = input[1] + c[0];
        input2[2] = input[2] + c[1];
        input2[3] = input[3] + c[2];

        std::vector<T> posOut(3,T());
        posOut = this->_superGeometry.getPhysR(globIC, input2[1], input2[2], input2[3]);
        _indicator(inside, &(posOut[0]) );
        if ( util::nearZero(inside[0]) ) {
          int overlap = this->_sLattice.getOverlap();
          // Get f_q of next fluid cell where l = opposite(q)
          T f = this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
          // Get f_l of the boundary cell
          // Add f_q and f_opp
          f += this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
          // Update force
          force[0] -= c[0]*f;
          force[1] -= c[1]*f;
          force[2] -= c[2]*f;
        }
      }
      force[0] = this->_converter.physForce(force[0]);
      force[1] = this->_converter.physForce(force[1]);
      force[2] = this->_converter.physForce(force[2]);

      output[0] = (posIn[1]-_indicator.getPos()[1])*force[2] - (posIn[2]-_indicator.getPos()[2])*force[1];
      output[1] = (posIn[2]-_indicator.getPos()[2])*force[0] - (posIn[0]-_indicator.getPos()[0])*force[2];
      output[2] = (posIn[0]-_indicator.getPos()[0])*force[1] - (posIn[1]-_indicator.getPos()[1])*force[0];

      return true;
    } else {
      return true;
    }
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  output[0] = 0.;
  output[1] = 0.;
  output[2] = 0.;
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    if (this->_superGeometry.get(input) == _material) {
      for (int iPop = 1; iPop < DESCRIPTOR<T>::q; ++iPop) {
        // Get direction
        const int* c = DESCRIPTOR<T>::c[iPop];
        // Get next cell located in the current direction
        // Check if the next cell is a fluid node
        if (this->_superGeometry.get(input[0], input[1] + c[0], input[2] + c[1],
                                     input[3] + c[2]) == 1) {
          int overlap = this->_sLattice.getOverlap();
          // Get f_q of next fluid cell where l = opposite(q)
          T f = this->_sLattice.getExtendedBlockLattice(
                  this->_sLattice.getLoadBalancer().loc(globIC)).get(
                  locix + overlap + c[0], lociy + overlap + c[1],
                  lociz + overlap + c[2])[iPop];
          // Get f_l of the boundary cell
          // Add f_q and f_opp
          f += this->_sLattice.getExtendedBlockLattice(
                 this->_sLattice.getLoadBalancer().loc(globIC)).get(
                 locix + overlap, lociy + overlap, lociz + overlap)[util::opposite<
                     DESCRIPTOR<T> >(iPop)];
          // Update force
          output[0] -= c[0] * (f - 2. * DESCRIPTOR<T>::t[iPop]);
          output[1] -= c[1] * (f - 2. * DESCRIPTOR<T>::t[iPop]);
          output[2] -= c[2] * (f - 2. * DESCRIPTOR<T>::t[iPop]);
        }
      }
      output[0] = this->_converter.physForce(output[0]);
      output[1] = this->_converter.physForce(output[1]);
      output[2] = this->_converter.physForce(output[2]);
      return true;
    } else {
      return true;
    }
  } else {
    return false;
  }
}

//end TODO

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeExternalField3D<T, DESCRIPTOR>::SuperLatticeExternalField3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, int beginsAt, int sizeOf)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, sizeOf), _beginsAt(beginsAt),
    _sizeOf(sizeOf)
{
  this->getName() = "externalField";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticeExternalField3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), _beginsAt, _sizeOf));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeExternalField3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePorosity3D<T, DESCRIPTOR>::SuperLatticePorosity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "porosity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockLatticePorosity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePorosity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPermeability3D<T, DESCRIPTOR>::SuperLatticePhysPermeability3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "permeability";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back( new BlockLatticePhysPermeability3D<T, DESCRIPTOR>(
                               this->_sLattice.getBlockLattice(iC), this->getConverter() ) );
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysPermeability3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::SuperLatticePhysCroppedPermeability3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "cropped_permeability";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); iC++ ) {
    this->_blockF.push_back( new BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>(
                               this->_sLattice.getBlockLattice(iC), this->getConverter() ) );
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysDarcyForce3D<T, DESCRIPTOR>::SuperLatticePhysDarcyForce3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDarcyForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  SuperLatticePhysPermeability3D<T, DESCRIPTOR> permeability(this->_sLattice, this->_converter);
  SuperLatticeVelocity3D<T, DESCRIPTOR> velocity(this->_sLattice);

  T nu = this->_converter.getCharNu();
  T K;
  T u[velocity.getTargetDim()];
  permeability(&K,input);
  velocity(u,input);

  output[0] = -nu / K * u[0];
  output[1] = -nu / K * u[1];
  output[2] = -nu / K * u[2];

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeAverage3D<T, DESCRIPTOR>::SuperLatticeAverage3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& superGeometry,
  const int material, T radius)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), f.getTargetDim()),
    _f(f),  _superGeometry(superGeometry), _material(material), _radius(radius)
{
  this->getName() = "Average(" + _f.getName() + ")";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeAverage3D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice().getLoadBalancer();

  //create boolean indicator functor isInSphere
  std::vector<T> center(3,T());
  cGeometry.getPhysR(&(center[0]), &(input[0]));
  IndicatorSphere3D<T> isInSphere(center, _radius);

  // iterate over all cuboids & points and test for material && isInSphere
  int numVoxels(0);
  if (this->_superGeometry.get(input) == _material) {
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      for (int iX = 0; iX < nX; ++iX) {
        for (int iY = 0; iY < nY; ++iY) {
          for (int iZ = 0; iZ < nZ; ++iZ) {
            std::vector<int> testLatticeR(input,input+4);
            //            testLatticeR[0] = load.glob(iC);
            //            testLatticeR[1] = iX;
            //            testLatticeR[2] = iY;
            //            testLatticeR[3] = iZ;
            T testPhysR[3] = {0};
            cGeometry.getPhysR(testPhysR,&(testLatticeR[0]));
            bool inside[1];
            isInSphere(inside,testPhysR);
            if (this->_superGeometry.get(input) == _material && inside[0]) {
              int inputTmp[4]= {load.glob(iC),iX,iY,iZ};
              T outputTmp[_f.getTargetDim()];
              _f(outputTmp,inputTmp);
              for ( int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
                output[iDim] += outputTmp[iDim];
              }
              numVoxels++;
            }
          }
        }
      }
    }

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
    for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
#ifdef PARALLEL_MODE_MPI
      singleton::mpi().reduceAndBcast(output[iDim], MPI_SUM);
#endif
      if (numVoxels > 0) {
        output[iDim] /= numVoxels;
      }
    }
  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperEuklidNorm3D<T, DESCRIPTOR>::SuperEuklidNorm3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 1), _f(f)
{
  this->getName() = "EuklidNorm(" + _f.getName() + ")";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.push_back(new BlockEuklidNorm3D<T, DESCRIPTOR>(f.getBlockF(iC)));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperEuklidNorm3D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::SuperLatticeInterpPhysVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& conv)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "InterpVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>* foo =
      new BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>(
      sLattice.getExtendedBlockLattice(iC),
      conv,
      &sLattice.getCuboidGeometry().get(this->_sLattice.getLoadBalancer().glob(iC)),
      sLattice.getOverlap());
    bLattices.push_back(foo);
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::~SuperLatticeInterpPhysVelocity3D()
{
  for ( typename std::vector<BlockLatticeInterpPhysVelocity3D<T,DESCRIPTOR>* >::iterator it = bLattices.begin() ; it != bLattices.end(); ++it) {
    delete (*it);
  }
  bLattices.clear();
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::operator()(T output[],
    const T input[], const int iC)
{
  bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(output, input);
}

//template<typename T, template<typename U> class DESCRIPTOR>
//bool SuperLatticeInterpVelocity3D<T, DESCRIPTOR>::operator() (T output[], const int input[])
//{
//  return true;
//}

//template<typename T, template<typename U> class DESCRIPTOR>
//SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::SuperLatticeInterpPhysVelocity3Degree3D(
//  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& conv)
//  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
//{
//  this->getName() = "Interp3DegreeVelocity";
//  int maxC = this->_sLattice.getLoadBalancer().size();
//  this->_blockF.reserve(maxC);
//  for (int iC = 0; iC < maxC; iC++) {
//    BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>* foo =
//      new BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>(
//      sLattice.getExtendedBlockLattice(iC),
//      conv,
//      &sLattice.getCuboidGeometry().get(this->_sLattice.getLoadBalancer().glob(iC)),
//      sLattice.getOverlap());
//    bLattices.push_back(foo);
//  }
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//void SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(T output[],
//    const T input[], const int iC)
//{
//  bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(output, input);
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::SuperLatticeInterpDensity3Degree3D(
//  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& conv)
//  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
//{
//  this->getName() = "Interp3DegreeDensity";
//  int maxC = this->_sLattice.getLoadBalancer().size();
//  this->_blockF.reserve(maxC);
//  for (int iC = 0; iC < maxC; iC++) {
//    BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>* foo =
//      new BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>(
//      sLattice.getExtendedBlockLattice(iC),
//      conv,
//      &sLattice.getCuboidGeometry().get(this->_sLattice.getLoadBalancer().glob(iC)),
//      sLattice.getOverlap());
//    bLattices.push_back(foo);
//  }
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//void SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(T output[],
//    const T input[], const int iC)
//{
//  bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(output, input);
//}

}  // end namespace olb

#endif
