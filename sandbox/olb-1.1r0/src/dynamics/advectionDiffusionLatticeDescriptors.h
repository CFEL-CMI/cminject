/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef ADVECTION_DIFFUSION_LATTICE_DESCRIPTORS_H
#define ADVECTION_DIFFUSION_LATTICE_DESCRIPTORS_H

#include <vector>
#include "latticeDescriptors.h"

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {
//===========================================================================//
//=================== AdvectionDiffusion Lattice Descriptors=================//
//===========================================================================//

/// D2Q5 lattice
template <typename T> struct D2Q5DescriptorBase {
  typedef D2Q5DescriptorBase<T> BaseDescriptor;
  enum { d = 2, q = 5 };      ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
};

struct Velocity2dDescriptor {
  static const int numScalars = 2;
  static const int numSpecies = 1;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 2;
};

struct Velocity2dBase {
  typedef Velocity2dDescriptor ExternalField;
};

struct Velocity3dDescriptor {
  static const int numScalars = 3;
  static const int numSpecies = 1;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 3;
};

struct Velocity3dBase {
  typedef Velocity3dDescriptor ExternalField;
};

struct particleAdvectionDiffusion3dDescriptor {
  static const int numScalars = 6;
  static const int numSpecies = 2;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 3;
  static const int velocity2BeginsAt = 3;
  static const int sizeOfVelocity2   = 3;
};

struct particleAdvectionDiffusion3dDescriptorBase {
  typedef particleAdvectionDiffusion3dDescriptor ExternalField;
};

/// AD D2Q5 lattice
template <typename T>
struct AdvectionDiffusionD2Q5Descriptor : public D2Q5DescriptorBase<T>, public Velocity2dBase {};

/// D3Q7 lattice
template <typename T>
struct D3Q7DescriptorBase {
  typedef D3Q7DescriptorBase<T> BaseDescriptor;
  enum { d = 3, q = 7 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound


};

template <typename T>
struct AdvectionDiffusionD3Q7Descriptor : public D3Q7DescriptorBase<T>, public Velocity3dBase {};

template <typename T>
struct particleAdvectionDiffusionD3Q7Descriptor : public D3Q7DescriptorBase<T>, public particleAdvectionDiffusion3dDescriptorBase { };


/// D3Q7 lattice for radiative transport problems @2016 A. Mink et al.
/// D3Q7 lattice
template <typename T>
struct D3Q7DescriptorBaseRTLB {
  typedef D3Q7DescriptorBaseRTLB<T> BaseDescriptor;
  enum { d = 3, q = 7 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const double henyeyPhaseFunction[q][q]; ///<anisotropic discrete scattering coefficient
};

// TODO: AM, Why inherit from Velocity3dBase?
template <typename T>
struct D3Q7DescriptorRTLB : public D3Q7DescriptorBaseRTLB<T>, public Velocity3dBase { };


/// D3Q19 lattice TODO: AM
template <typename T>
struct D3Q19DescriptorBaseRTLB {
  typedef D3Q19DescriptorBaseRTLB<T> BaseDescriptor;
  enum { d = 3, q = 19 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const double henyeyPhaseFunction[q][q]; ///<anisotropic discrete scattering coefficient
};

template <typename T>
struct D3Q19DescriptorRTLB  : public D3Q19DescriptorBaseRTLB<T>, public Velocity3dBase {};



/// D3Q27 lattice
template <typename T>
struct D3Q27DescriptorBaseRTLB {
  typedef D3Q27DescriptorBaseRTLB<T> BaseDescriptor;
  enum { d = 3, q = 27 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const double henyeyPhaseFunction[q][q]; ///<anisotropic discrete scattering coefficient
};


template <typename T>
struct D3Q27DescriptorRTLB  : public D3Q27DescriptorBaseRTLB<T>, public Velocity3dBase {};


}  // namespace descriptors

}  // namespace olb

#endif
