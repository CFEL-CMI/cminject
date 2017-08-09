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

#ifndef BLOCK_CALC_F_3D_H
#define BLOCK_CALC_F_3D_H


#include "functors/blockBaseF3D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// arithmetic helper class for BlockF3D functors
template <typename T>
class BlockCalc3D : public BlockF3D<T> {
protected:
  BlockF3D<T>& _f;
  BlockF3D<T>& _g;
public:
  BlockCalc3D(BlockF3D<T>& f, BlockF3D<T>& g);
};

/// addition functor
template <typename T>
class BlockPlus3D : public BlockCalc3D<T> {
public:
  BlockPlus3D(BlockF3D<T>& f, BlockF3D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// subtraction functor
template <typename T>
class BlockMinus3D : public BlockCalc3D<T> {
public:
  BlockMinus3D(BlockF3D<T>& f, BlockF3D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// multiplication functor
template <typename T>
class BlockMultiplication3D : public BlockCalc3D<T> {
public:
  BlockMultiplication3D(BlockF3D<T>& f, BlockF3D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// division functor
template <typename T>
class BlockDivision3D : public BlockCalc3D<T> {
public:
  BlockDivision3D(BlockF3D<T>& f,BlockF3D<T>& g);
  bool operator() (T output[], const int input[]);
};


} // end namespace olb

#endif
