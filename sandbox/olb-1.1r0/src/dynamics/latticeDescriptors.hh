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
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- generic code
 */
#ifndef LATTICE_DESCRIPTORS_HH
#define LATTICE_DESCRIPTORS_HH

#include "latticeDescriptors.h"

namespace olb {

namespace descriptors {

// D2Q9 ////////////////////////////////////////////////////////////

template<typename T>
const int D2Q9DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D2Q9DescriptorBase<T>::c
[D2Q9DescriptorBase<T>::q][D2Q9DescriptorBase<T>::d] = {
  { 0, 0},
  {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
  { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
};

template<typename T>
const int D2Q9DescriptorBase<T>::opposite[D2Q9DescriptorBase<T>::q] = {
  0, 5, 6, 7, 8, 1, 2, 3, 4
};


template<typename T>
const T D2Q9DescriptorBase<T>::t[D2Q9DescriptorBase<T>::q] = {
  (T)4/(T)9, (T)1/(T)36, (T)1/(T)9, (T)1/(T)36, (T)1/(T)9,
  (T)1/(T)36, (T)1/(T)9, (T)1/(T)36, (T)1/(T)9
};

template<typename T>
const T D2Q9DescriptorBase<T>::invCs2 = (T)3;

// D3Q13 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q13DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q13DescriptorBase<T>::c
[D3Q13DescriptorBase<T>::q][D3Q13DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template<typename T>
const int D3Q13DescriptorBase<T>::opposite[D3Q13DescriptorBase<T>::q] = {
  0, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6
};

template<typename T>
const T D3Q13DescriptorBase<T>::t[D3Q13DescriptorBase<T>::q] = {
  (T)1/(T)2,

  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,
  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,

  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,
  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24
};


/** This parameter is chosen to enhance numerical stability */
template<typename T>
const T D3Q13DescriptorBase<T>::invCs2 = (T)3;

/** This parameter is chosen to enhance numerical stability */
template<typename T>
const T D3Q13DescriptorBase<T>::lambda_e = (T)1.5;

/** This parameter is chosen to enhance numerical stability */
template<typename T>
const T D3Q13DescriptorBase<T>::lambda_h = (T)1.8;


// D3Q15 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q15DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q15DescriptorBase<T>::c
[D3Q15DescriptorBase<T>::q][D3Q15DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1},

};

template<typename T>
const int D3Q15DescriptorBase<T>::opposite[D3Q15DescriptorBase<T>::q] = {
  0, 8, 9, 10, 11, 12, 13, 14, 1, 2, 3, 4, 5, 6, 7
};

template<typename T>
const T D3Q15DescriptorBase<T>::t[D3Q15DescriptorBase<T>::q] = {
  (T)2/(T)9,

  (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
  (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72,

  (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
  (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72
};

template<typename T>
const T D3Q15DescriptorBase<T>::invCs2 = (T)3;


// D3Q19 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q19DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q19DescriptorBase<T>::c
[D3Q19DescriptorBase<T>::q][D3Q19DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template<typename T>
const int D3Q19DescriptorBase<T>::opposite[D3Q19DescriptorBase<T>::q] = {
  0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9
};

template<typename T>
const T D3Q19DescriptorBase<T>::t[D3Q19DescriptorBase<T>::q] = {
  (T)1/(T)3,

  (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,

  (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36
};

template<typename T>
const T D3Q19DescriptorBase<T>::invCs2 = (T)3;


// D3Q27 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q27DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q27DescriptorBase<T>::c
[D3Q27DescriptorBase<T>::q][D3Q27DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
};

template<typename T>
const int D3Q27DescriptorBase<T>::opposite[D3Q27DescriptorBase<T>::q] = {
  0, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
};

template<typename T>
const T D3Q27DescriptorBase<T>::t[D3Q27DescriptorBase<T>::q] = {
  (T)8/(T)27,

  (T)2/(T)27, (T)2/(T)27, (T)2/(T)27,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)216, (T)1/(T)216, (T)1/(T)216, (T)1/(T)216,

  (T)2/(T)27, (T)2/(T)27, (T)2/(T)27,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)216, (T)1/(T)216, (T)1/(T)216, (T)1/(T)216
};

template<typename T>
const T D3Q27DescriptorBase<T>::invCs2 = (T)3;

}  // namespace descriptors

}  // namespace olb

#endif
