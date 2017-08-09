/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathen
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

#ifndef SHEAR_SMAGORINSKY_LATTICE_DESCRIPTOR_H
#define SHEAR_SMAGORINSKY_LATTICE_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {

// 2D Descriptors for flow with Shear-Improved Smagorinsky

struct ShearSmagorinsky2dDescriptor {
  static const int numScalars = 1;
  static const int numSpecies = 1;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
};

struct ShearSmagorinsky2dDescriptorBase {
  typedef ShearSmagorinsky2dDescriptor ExternalField;
};

template <typename T> struct ShearSmagorinskyD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ShearSmagorinsky2dDescriptorBase {
};


/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Smagorinsky

struct ShearSmagorinsky3dDescriptor {
  static const int numScalars = 1;
  static const int numSpecies = 1;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
};

struct ShearSmagorinsky3dDescriptorBase {
  typedef ShearSmagorinsky3dDescriptor ExternalField;
};

template <typename T> struct ShearSmagorinskyD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ShearSmagorinsky3dDescriptorBase {
};


/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Forced Shear-Improved Smagorinsky

struct ForcedShearSmagorinsky3dDescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 2;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
  static const int forceBeginsAt    = 1;
  static const int sizeOfForce      = 3;
};

struct ForcedShearSmagorinsky3dDescriptorBase {
  typedef ForcedShearSmagorinsky3dDescriptor ExternalField;
};

template <typename T> struct ForcedShearSmagorinskyD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ForcedShearSmagorinsky3dDescriptorBase {
};

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Forced Shear-Improved Smagorinsky

struct ForcedShearWallSmagorinsky3dDescriptor {
  static const int numScalars = 5;
  static const int numSpecies = 3;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
  static const int forceBeginsAt    = 1;
  static const int sizeOfForce      = 3;
  static const int tauWIsAt   = 4;
  static const int sizeOfTauW      = 1;
};

struct ForcedShearWallSmagorinsky3dDescriptorBase {
  typedef ForcedShearWallSmagorinsky3dDescriptor ExternalField;
};

template <typename T> struct ForcedShearWallSmagorinskyD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ForcedShearWallSmagorinsky3dDescriptorBase {
};


} // namespace descriptors

} // namespace olb

#endif
