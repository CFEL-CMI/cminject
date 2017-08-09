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
#ifndef BLOCK_CALC_F_3D_HH
#define BLOCK_CALC_F_3D_HH


#include "functors/blockCalcF3D.h"


namespace olb {


template <typename T>
BlockCalc3D<T>::BlockCalc3D (BlockF3D<T>& f, BlockF3D<T>& g)
  : BlockF3D<T>( f.getBlockStructure(), f.getTargetDim() ), _f(f), _g(g) { }


template <typename T>
BlockPlus3D<T>::BlockPlus3D(BlockF3D<T>& f, BlockF3D<T>& g)
  : BlockCalc3D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockPlus3D<T>::operator()(T output[], const int input[])
{
  this->_g(output,input);
  T tmp[this->_f.getTargetDim()];
  this->_f(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] += tmp[i];
  }
  return true;
}


template <typename T>
BlockMinus3D<T>::BlockMinus3D(BlockF3D<T>& f, BlockF3D<T>& g)
  : BlockCalc3D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockMinus3D<T>::operator()(T output[], const int input[])
{
  this->_g(output,input);
  T tmp[this->_f.getTargetDim()];
  this->_f(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] -= tmp[i];
  }
  return true;
}


template <typename T>
BlockMultiplication3D<T>::BlockMultiplication3D(BlockF3D<T>& f, BlockF3D<T>& g)
  : BlockCalc3D<T>(f,g)
{
  this->getName() = f.getName() + "*" + g.getName();
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockMultiplication3D<T>::operator()(T output[], const int input[])
{
  this->_g(output,input);
  T tmp[this->_f.getTargetDim()];
  this->_f(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] *= tmp[i];
  }
  return true;
}


template <typename T>
BlockDivision3D<T>::BlockDivision3D (BlockF3D<T>& f, BlockF3D<T>& g)
  : BlockCalc3D<T>(f,g)
{
  this->getName() = f.getName() + "/" + g.getName();
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockDivision3D<T>::operator()(T output[], const int input[])
{
  this->_g(output,input);
  T tmp[this->_f.getTargetDim()];
  this->_f(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] /= tmp[i];
  }
  return true;
}



/////////////////////////////////operator()/// ////////////////////////////////
template <typename T>
BlockF3D<T>& BlockF3D<T>::operator+(BlockF3D<T>& rhs)
{
  auto tmp = std::make_shared< BlockPlus3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF3D<T>& BlockF3D<T>::operator-(BlockF3D<T>& rhs)
{
  auto tmp = std::make_shared< BlockMinus3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF3D<T>& BlockF3D<T>::operator*(BlockF3D<T>& rhs)
{
  auto tmp = std::make_shared< BlockMultiplication3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF3D<T>& BlockF3D<T>::operator/(BlockF3D<T>& rhs)
{
  auto tmp = std::make_shared< BlockDivision3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


/////////////////// BlockDataF3D ///////////////////////////

template <typename T,typename BaseType>
BlockDataF3D<T,BaseType>::BlockDataF3D(BlockData3D<T,BaseType>& blockData)
  : BlockF3D<T>( blockData, blockData.getSize() ), _blockData(blockData),
    _isConstructed(false)
{}

template <typename T,typename BaseType>
BlockDataF3D<T,BaseType>::BlockDataF3D(BlockF3D<BaseType>& f)
  : BlockF3D<T>( f.getBlockStructure(), f.getTargetDim() ),
    _blockData( *(new BlockData3D<T,BaseType>(f)) ), _isConstructed(true)
{}


template <typename T,typename BaseType>
BlockDataF3D<T,BaseType>::BlockDataF3D(int nx, int ny, int nz, int size)
  : BlockF3D<T>( *(new BlockData3D<T,BaseType>(nx, ny, nz, size)), size ), _blockData( static_cast<BlockData3D<T,BaseType>&>(this->getBlockStructure() )),
    _isConstructed(true)
{}

template <typename T,typename BaseType>
BlockDataF3D<T,BaseType>::~BlockDataF3D()
{
  if (_isConstructed) {
    delete &_blockData;
  }
}

template <typename T,typename BaseType>
BlockData3D<T,BaseType>& BlockDataF3D<T,BaseType>::getBlockData()
{
  return _blockData;
}

// access _blockLattice3D
template <typename T, typename BaseType>
bool BlockDataF3D<T, BaseType>::operator()(BaseType output[], const int input[])
{
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = _blockData.get( input[0], input[1], input[2], iDim );
  }
  return true;
}



} // end namespace olb

#endif
