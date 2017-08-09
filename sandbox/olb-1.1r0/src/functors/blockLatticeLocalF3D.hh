/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_LOCAL_F_3D_HH
#define BLOCK_LATTICE_LOCAL_F_3D_HH


#include<cmath>

#include "functors/blockLatticeLocalF3D.h"
#include "functors/blockBaseF3D.h"
#include "core/blockLatticeStructure3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeFpop3D<T, DESCRIPTOR>::BlockLatticeFpop3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR<T>::q)
{
  this->getName() = "fPop";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeFpop3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int iPop = 0; iPop < DESCRIPTOR<T>::q; ++iPop) {
    output[iPop] =
      this->_blockLattice.get(input[0], input[1], input[2])[iPop]
      + DESCRIPTOR<T>::t[iPop];
  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeDissipation3D<T, DESCRIPTOR>::BlockLatticeDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter)
{
  this->getName() = "dissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = _converter.getLatticeNu();
  T omega = _converter.getOmega();
  output[0] = PiNeqNormSqr * nuLattice
              * pow(omega * DESCRIPTOR<T>::invCs2, 2) / rho / 2.;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysDissipation3D<T, DESCRIPTOR>::BlockLatticePhysDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter)
{
  this->getName() = "physDissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = _converter.getLatticeNu();
  T omega = _converter.getOmega();
  T dt = _converter.physTime();
  output[0] = PiNeqNormSqr * nuLattice
              * pow(omega * DESCRIPTOR<T>::invCs2 / rho, 2) / 2.
              * _converter.getCharNu() / _converter.getLatticeNu() / dt / dt;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticeEffevtiveDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter, T smagoConst)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst)
{
  this->getName() = "EffevtiveDissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = _converter.getLatticeNu();
  T omega = _converter.getOmega();

  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega;
  /// Turbulent realaxation time
  //T dx = _converter.getLatticeL();
  //T dt = _converter.getDeltaT();
  T preFactor = (_smagoConst*_smagoConst)*DESCRIPTOR<T>::invCs2*DESCRIPTOR<T>::invCs2*2*sqrt(2);
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor/rho*PiNeqNorm) - tau_mol);
  //T omegaTurbulent = 1./tau_turb;
  T omegaEff = 1./(tau_turb + tau_mol);
  T nuTurbulentLattice = tau_turb*DESCRIPTOR<T>::invCs2;

  output[0] = PiNeqNormSqr * (nuTurbulentLattice + nuLattice)
              * pow(omegaEff * DESCRIPTOR<T>::invCs2, 2) / rho / 2.;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticePhysEffevtiveDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter, T smagoConst)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst)
{
  this->getName() = "physEffevtiveDissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = _converter.getLatticeNu();
  T omega = _converter.getOmega();

  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega;
  /// Turbulent realaxation time
  T dt = _converter.getDeltaT();
  T preFactor = (_smagoConst*_smagoConst)*DESCRIPTOR<T>::invCs2*DESCRIPTOR<T>::invCs2*2*sqrt(2);
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor/rho*PiNeqNorm) - tau_mol);
//  T omegaTurbulent = 1./tau_turb;
  T omegaEff = 1./(tau_turb + tau_mol);
  T nuTurbulentLattice = tau_turb/DESCRIPTOR<T>::invCs2;
  //std::cout << nuTurbulentLattice << std::endl;

  output[0] = PiNeqNormSqr * (nuTurbulentLattice + nuLattice)
              * pow(omegaEff * DESCRIPTOR<T>::invCs2 / rho, 2) / 2.
              * _converter.getCharNu() / _converter.getLatticeNu() / dt / dt;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeDensity3D<T, DESCRIPTOR>::BlockLatticeDensity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "density";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeDensity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = this->_blockLattice.get(input[0], input[1], input[2]).computeRho();
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeVelocity3D<T, DESCRIPTOR>::BlockLatticeVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3)
{
  this->getName() = "velocity";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho;
  this->_blockLattice.get(input[0], input[1], input[2]).computeRhoU(rho, output);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeStrainRate3D<T, DESCRIPTOR>::BlockLatticeStrainRate3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9)
{
  this->getName() = "strainRate";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T omega = this->_converter.getOmega();

  output[0] = -pi[0] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[1] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[2] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[3] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[4] = -pi[3] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[5] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[6] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[7] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[8] = -pi[5] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::BlockLatticePhysStrainRate3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9)
{
  this->getName() = "strainRate";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T omega = this->_converter.getOmega();
  T dt = this->_converter.physTime();

  output[0] = -pi[0] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[1] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[2] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[3] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[4] = -pi[3] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[5] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[6] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[7] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[8] = -pi[5] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeGeometry3D<T, DESCRIPTOR>::BlockLatticeGeometry3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
    BlockGeometryStructure3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "geometry";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeGeometry3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = _blockGeometry.getMaterial(input[0], input[1], input[2]);

  if (_material != -1) {
    if ( util::nearZero(_material-output[0]) ) {
      output[0] = 1.;
      return true;
    } else {
      output[0] = 0.;
      return true;
    }
  }
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeRank3D<T,DESCRIPTOR>::BlockLatticeRank3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "rank";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeRank3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = singleton::mpi().getRank() + 1;
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeCuboid3D<T,DESCRIPTOR>::BlockLatticeCuboid3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int iC)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _iC(iC)
{
  this->getName() = "cuboid";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeCuboid3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = _iC + 1;
  return true;
}



template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysPressure3D<T, DESCRIPTOR>::BlockLatticePhysPressure3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "physPressure";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysPressure3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho = this->_blockLattice.get(input[0], input[1], input[2]).computeRho();
  output[0]=this->_converter.physPressureFromRho(rho);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysVelocity3D<T, DESCRIPTOR>::BlockLatticePhysVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter, bool print)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _print(print)
{
  this->getName() = "physVelocity";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho;
  if (_print) {
    std::cout << input[0] << " " << input[1] << " " << input[2] << " | "
              << singleton::mpi().getRank() << std::endl;
  }
  this->_blockLattice.get(input[0], input[1], input[2]).computeRhoU(rho, output);
  output[0] = this->_converter.physVelocity(output[0]);
  output[1] = this->_converter.physVelocity(output[1]);
  output[2] = this->_converter.physVelocity(output[2]);

  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalPorosity3D<T,DESCRIPTOR>::BlockLatticePhysExternalPorosity3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtPorosityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalPorosity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(input[0],input[1],input[2]).computeExternalField(DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, output);
//  output[0] > 0 ? cout << output[0] << std::endl : cout; ;
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalVelocity3D<T,DESCRIPTOR>::BlockLatticePhysExternalVelocity3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtVelocityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalVelocity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(input[0],input[1],input[2]).computeExternalField(1, 2, output);
  //  this->_blockLattice.get(input[0],input[1]).computeExternalField(DESCRIPTOR<T>::ExternalField::localDragBeginsAt,
  //                                                                  DESCRIPTOR<T>::ExternalField::sizeOfLocalDrag, output);
  output[0]=this->_converter.physVelocity(output[0]);
  output[1]=this->_converter.physVelocity(output[1]);
  output[2]=this->_converter.physVelocity(output[2]);

  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::BlockLatticePhysExternalParticleVelocity3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtParticleVelocityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T* foo = this->_blockLattice.get(input[0],input[1],input[2]).getExternal(0);

  if (foo[4] > std::numeric_limits<T>::epsilon()) {
    output[0]=this->_converter.physVelocity(foo[1]/foo[4]);
    output[1]=this->_converter.physVelocity(foo[2]/foo[4]);
    output[2]=this->_converter.physVelocity(foo[3]/foo[4]);
    return true;
  }
  output[0]=this->_converter.physVelocity(foo[1]);
  output[1]=this->_converter.physVelocity(foo[2]);
  output[2]=this->_converter.physVelocity(foo[3]);
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeExternal3D<T, DESCRIPTOR>::BlockLatticeExternal3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, int start, int size)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, size), _start(start), _size(size)
{
  this->getName() = "extField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeExternal3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    _start, _size, output);
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysExternal3D<T, DESCRIPTOR>::BlockLatticePhysExternal3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3)
{
  this->getName() = "physExtField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysExternal3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::velocityBeginsAt,
    DESCRIPTOR<T>::ExternalField::sizeOfVelocity, output);
  output[0] = this->_converter.physVelocity(output[0]);
  output[1] = this->_converter.physVelocity(output[1]);
  output[2] = this->_converter.physVelocity(output[2]);
  return true;
}



template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physBoundaryForce";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysBoundaryForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();
  output[1] = T();
  output[2] = T();

  if (this->_blockGeometry.get(input[0],input[1],input[2]) == _material) {
    for (int iPop = 1; iPop < DESCRIPTOR<T>::q; ++iPop) {
      // Get direction
      const int* c = DESCRIPTOR<T>::c[iPop];
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get(input[0] + c[0], input[1] + c[1], input[2] + c[2]) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + c[0], input[1] + c[1], input[2] + c[2])[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input[0], input[1], input[2])[util::opposite<DESCRIPTOR<T> >(iPop)];
        // Update force
        output[0] -= c[0] * f;
        output[1] -= c[1] * f;
        output[2] -= c[2] * f;
      }
    }
    output[0] = this->_converter.physForce(output[0]);
    output[1] = this->_converter.physForce(output[1]);
    output[2] = this->_converter.physForce(output[2]);
    return true;
  } else {
    return true;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR<T>::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (BlockGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->_blockLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[1] -= c[1]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[2] -= c[2]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //        }
  //      }
  //      force[0]=this->_converter.physForce(force[0]);
  //      force[1]=this->_converter.physForce(force[1]);
  //      force[2]=this->_converter.physForce(force[2]);
  //      return force;
  //    }
  //    else {
  //      return force;
  //    }
  //  }
  //  else {
  //    return force;
  //  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeExternalField3D<T, DESCRIPTOR>::BlockLatticeExternalField3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, int beginsAt, int sizeOf)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, sizeOf),
    _beginsAt(beginsAt), _sizeOf(sizeOf)
{
  this->getName() = "externalField";
}

// under construction
template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeExternalField3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(_beginsAt, _sizeOf, output);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePorosity3D<T, DESCRIPTOR>::BlockLatticePorosity3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "porosity";
}

// under construction
template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePorosity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::porosityIsAt, DESCRIPTOR<T>::ExternalField::sizeOfPorosity, output);
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysPermeability3D<T, DESCRIPTOR>::BlockLatticePhysPermeability3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter)
{
  this->getName() = "permeability";
}

//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticePhysPermeability3D<T, DESCRIPTOR>::BlockLatticePhysPermeability3D(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
//  BlockGeometry3D<T>& blockGeometry, int material,
//  const LBconverter<T>& converter)
//  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1),
//    _blockGeometry(blockGeometry),
//    _material(material)
//{
//  this->getName() = "permeability";
//}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysPermeability3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T value;
  // get porosity
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::porosityIsAt, DESCRIPTOR<T>::ExternalField::sizeOfPorosity, &value);
  // convert to physPermeability
  output[0]=_converter.physPermeability(value);
  // \todo Why do we need this???
  if (output[0] >= 42 && output[0] <= 42 && output[0] != 42) {
    output[0] = 999999;
  }
  if (std::isinf(output[0])) {
    output[0] = 1e6;
  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::BlockLatticePhysCroppedPermeability3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter)
{
  this->getName() = "cropped_permeability";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T value;
  // get porosity
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::porosityIsAt, DESCRIPTOR<T>::ExternalField::sizeOfPorosity, &value);
  // convert to physPermeability
  output[0]=_converter.physPermeability(value);
  // \todo Why do we need this???
  if (output[0] >= 42 && output[0] <= 42 && output[0] != 42) {
    output[0] = 1;
  }
  if (std::isinf(output[0])) {
    output[0] = 1;
  }
  if (output[0] > 1.) {
    output[0] = 1.;
  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::BlockLatticePhysDarcyForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "alphaU";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  BlockLatticePhysPermeability3D<T, DESCRIPTOR> permeability(this->_blockLattice, this->_converter);
  BlockLatticeVelocity3D<T, DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getCharNu();
  permeability(output,input);
  T K = output[0];
  velocity(output,input);

  output[0] *= -nu / K;
  output[1] *= -nu / K;
  output[2] *= -nu / K;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeAverage3D<T, DESCRIPTOR>::BlockLatticeAverage3D(BlockLatticeF3D<T, DESCRIPTOR>& f,
    BlockGeometry3D<T>& blockGeometry, int material, T radius)
  : BlockLatticeF3D<T, DESCRIPTOR>(f.getBlockLattice(), f.getTargetDim()),
    _f(f),
    _blockGeometry(blockGeometry),
    _material(material),
    _radius(radius)
{
  this->getName() = "Average(" + f.getName() + ")";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeAverage3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  //  CuboidGeometry3D<T>& cGeometry = f.getBlockLattice().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice().get_load();

  //  //create boolean indicator functor isInSphere
  //  std::vector<T> center(3,0);
  //  center[0] = (int)cGeometry.get(load.glob(input[0])).get_globPosX() + input[1];
  //  center[1] = (int)cGeometry.get(load.glob(input[0])).get_globPosY() + input[2];
  //  center[2] = (int)cGeometry.get(load.glob(input[0])).get_globPosZ() + input[3];
  //  SphereAnalyticalF3D<bool,T> isInSphere(center,radius);

  // iterate over all cuboids & points and test for material && isInSphere
  output[0]=0;
  //  int numVoxels(0);
  //  if (this->blockGeometry.getMaterial(center[0],center[1],center[2]) == material) {
  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      for (int iX=0; iX<nX; ++iX) {
  //        for (int iY=0; iY<nY; ++iY) {
  //          for (int iZ=0; iZ<nZ; ++iZ) {
  //            std::vector<T> glob(3,0);
  //            glob[0] = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            glob[1] = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            glob[2] = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->blockGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
  //                && isInSphere(glob)[0]==true) {
  //              for (unsigned iD=0; iD<f(load.glob(0),0,0,0).size(); iD++) {
  //                tmp[iD]+=f(load.glob(iC),iX,iY,iZ)[iD];
  //              }
  //              numVoxels++;
  //            }
  //          }
  //        }
  //      }
  //    }

  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  //#endif
  //    for (int iD=0; iD<f.getTargetDim(); iD++) {
  //#ifdef PARALLEL_MODE_MPI
  //      singleton::mpi().reduceAndBcast(tmp[iD], MPI_SUM);
  //#endif
  //      if (numVoxels>0) {
  //        tmp[iD] /= numVoxels;
  //      }
  //    }
  //  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockEuklidNorm3D<T, DESCRIPTOR>::BlockEuklidNorm3D(BlockF3D<T>& f)
  : BlockF3D<T>(f.getBlockStructure(), 1),
    _f(f)
{
  this->getName() = "EuklidNorm(" + f.getName() + ")";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockEuklidNorm3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();  // flash output, since this methods adds values.
  T data[_f.getTargetDim()];
  _f(data,input);
  for (int i = 0; i < _f.getTargetDim(); ++i) {
    output[0] += data[i] * data[i];
  }
  output[0] = sqrt(output[0]);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, LBconverter<T>& conv, Cuboid3D<T>* c, int overlap)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
    _conv(conv),
    _cuboid(c),
    _overlap(overlap)
{
  this->getName() = "BlockLatticeInterpVelocity3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3D(
  const BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& rhs) :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
  _conv(rhs._conv),
  _cuboid(rhs._cuboid)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
void BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::operator()(T output[3], const T input[3])
{
  T u[3], rho, volume;
  T d[3], e[3];
  int latIntPos[3] = {0};
  T latPhysPos[3] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

  T deltaRinv = 1. / _cuboid->getDeltaR();
  d[0] = (input[0] - latPhysPos[0]) * deltaRinv;
  d[1] = (input[1] - latPhysPos[1]) * deltaRinv;
  d[2] = (input[2] - latPhysPos[2]) * deltaRinv;

  e[0] = 1. - d[0];
  e[1] = 1. - d[1];
  e[2] = 1. - d[2];

  latIntPos[0]+=_overlap;
  latIntPos[1]+=_overlap;
  latIntPos[2]+=_overlap;

  this->_blockLattice.get(latIntPos[0], latIntPos[1],
                          latIntPos[2]).computeRhoU(rho, u);
  volume = e[0] * e[1] * e[2];
  output[0] = u[0] * volume;
  output[1] = u[1] * volume;
  output[2] = u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1] + 1,
                          latIntPos[2]).computeRhoU(rho, u);
  volume = e[0] * d[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1],
                          latIntPos[2]).computeRhoU(rho, u);
  volume = d[0] * e[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1] + 1,
                          latIntPos[2]).computeRhoU(rho, u);
  volume = d[0] * d[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1],
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = e[0] * e[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1] + 1,
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = e[0] * d[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1],
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = d[0] * e[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1] + 1,
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = d[0] * d[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  output[0] = _conv.physVelocity(output[0]);
  output[1] = _conv.physVelocity(output[1]);
  output[2] = _conv.physVelocity(output[2]);
}

//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3Degree3D(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, LBconverter<T>& conv, Cuboid3D<T>* c, int overlap)
//  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
//    _conv(conv),
//    _cuboid(c),
//    _overlap(overlap)
//{
//  this->getName() = "BlockLatticeInterpVelocity3Degree3D";
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3Degree3D(
//  const BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>& rhs) :
//  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
//  _conv(rhs._conv),
//  _cuboid(rhs._cuboid)
//{
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//void BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(T output[3], const T input[3])
//{
//  T u[3], rho, volume;
//  int latIntPos[3] = {0};
//  T latPhysPos[3] = {T()};
//  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
//  _cuboid->getPhysR(latPhysPos, latIntPos);
//
//  latIntPos[0]+=_overlap;
//  latIntPos[1] += _overlap;
//  latIntPos[2] += _overlap;
//
//  volume=T(1);
//  for (int i = -1; i <= 2; ++i) {
//    for (int j = -1; j <= 2; ++j) {
//      for (int k = -1; k <= 2; ++k) {
//
//        this->_blockLattice.get(latIntPos[0]+i, latIntPos[1]+j, latIntPos[2]+k).computeRhoU(
//          rho, u);
//        for (int l = -1; l <= 2; ++l) {
//          if (l != i) {
//            volume *= (input[0] - (latPhysPos[0]+ l *_cuboid->getDeltaR()))
//                      / (latPhysPos[0] + i *_cuboid->getDeltaR()
//                         - (latPhysPos[0] + l *_cuboid->getDeltaR()));
//          }
//        }
//        for (int m = -1; m <= 2; ++m) {
//          if (m != j) {
//            volume *= (input[1]
//                       - (latPhysPos[1] + m *_cuboid->getDeltaR()))
//                      / (latPhysPos[1] + j * _cuboid->getDeltaR()
//                         - (latPhysPos[1] + m * _cuboid->getDeltaR()));
//          }
//        }
//        for (int n = -1; n <= 2; ++n) {
//          if (n != k) {
//            volume *= (input[2]
//                       - (latPhysPos[2] + n * _cuboid->getDeltaR()))
//                      / (latPhysPos[2] + k * _cuboid->getDeltaR()
//                         - (latPhysPos[2] + n * _cuboid->getDeltaR()));
//          }
//        }
//        output[0] += u[0] * volume;
//        output[1] += u[1] * volume;
//        output[2] += u[2] * volume;
//        volume=T(1);
//      }
//    }
//  }
//
//  output[0] = _conv.physVelocity(output[0]);
//  output[1] = _conv.physVelocity(output[1]);
//  output[2] = _conv.physVelocity(output[2]);
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, LBconverter<T>& conv, Cuboid3D<T>* c, int overlap)
//  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
//    _conv(conv),
//    _cuboid(c),
//    _overlap(overlap)
//{
//  this->getName() = "BlockLatticeInterpDensity3Degree3D";
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
//  const BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>& rhs) :
//  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
//  _conv(rhs._conv),
//  _cuboid(rhs._cuboid)
//{
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//void BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(T output[DESCRIPTOR<T>::q], const T input[3])
//{
//  T u[3], rho, volume, density;
//  int latIntPos[3] = {0};
//  T latPhysPos[3] = {T()};
//  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
//  _cuboid->getPhysR(latPhysPos, latIntPos);
//
//  latIntPos[0] += _overlap;
//  latIntPos[1] += _overlap;
//  latIntPos[2] += _overlap;
//
//  for (unsigned iPop = 0; iPop < DESCRIPTOR<T>::q; ++iPop) {
//    volume=T(1);
//    for (int i = -1; i <= 2; ++i) {
//      for (int j = -1; j <= 2; ++j) {
//        for (int k = -1; k <= 2; ++k) {
//
//          density = this->_blockLattice.get(latIntPos[0]+i, latIntPos[1]+j, latIntPos[2]+k).operator[](iPop);
//          for (int l = -1; l <= 2; ++l) {
//            if (l != i) {
//              volume *= (input[0] - (latPhysPos[0]+ l *_cuboid->getDeltaR()))
//                        / (latPhysPos[0] + i *_cuboid->getDeltaR()
//                           - (latPhysPos[0] + l *_cuboid->getDeltaR()));
//            }
//          }
//          for (int m = -1; m <= 2; ++m) {
//            if (m != j) {
//              volume *= (input[1]
//                         - (latPhysPos[1] + m *_cuboid->getDeltaR()))
//                        / (latPhysPos[1] + j * _cuboid->getDeltaR()
//                           - (latPhysPos[1] + m * _cuboid->getDeltaR()));
//            }
//          }
//          for (int n = -1; n <= 2; ++n) {
//            if (n != k) {
//              volume *= (input[2]
//                         - (latPhysPos[2] + n * _cuboid->getDeltaR()))
//                        / (latPhysPos[2] + k * _cuboid->getDeltaR()
//                           - (latPhysPos[2] + n * _cuboid->getDeltaR()));
//            }
//          }
//          output[iPop] += density * volume;
//
//          volume=T(1);
//        }
//      }
//    }
//  }
//}

}  // end namespace olb

#endif
