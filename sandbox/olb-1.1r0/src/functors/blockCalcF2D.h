/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Albert Mink, Lukas Baron, Mathias J. Krause
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

#ifndef BLOCK_CALC_F_2D_H
#define BLOCK_CALC_F_2D_H


#include "functors/blockBaseF2D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// arithmetic helper class for BlockLatticeF2D functors
template <typename T>
class BlockCalc2D : public BlockF2D<T> {
protected:
  BlockF2D<T>& _f;
  BlockF2D<T>& _g;
public:
  BlockCalc2D(BlockF2D<T>& f, BlockF2D<T>& g);
};

/// addition functor
template <typename T>
class BlockPlus2D : public BlockCalc2D<T> {
public:
  BlockPlus2D(BlockF2D<T>& f, BlockF2D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// subtraction functor
template <typename T>
class BlockMinus2D : public BlockCalc2D<T> {
public:
  BlockMinus2D(BlockF2D<T>& f, BlockF2D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// multiplication functor
template <typename T>
class BlockMultiplication2D : public BlockCalc2D<T> {
public:
  BlockMultiplication2D(BlockF2D<T>& f, BlockF2D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// division functor
template <typename T>
class BlockDivision2D : public BlockCalc2D<T> {
public:
  BlockDivision2D(BlockF2D<T>& f, BlockF2D<T>& g);
  bool operator() (T output[], const int input[]);
};


} // end namespace olb

#endif
