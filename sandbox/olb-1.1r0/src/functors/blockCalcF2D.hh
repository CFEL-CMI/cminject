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
#ifndef BLOCK_CALC_F_2D_HH
#define BLOCK_CALC_F_2D_HH


#include "functors/blockCalcF2D.h"


namespace olb {


template <typename T>
BlockCalc2D<T>::BlockCalc2D (BlockF2D<T>& f, BlockF2D<T>& g)
  : BlockF2D<T>( f.getBlockStructure(), f.getTargetDim() ), _f(f), _g(g) { }


template <typename T>
BlockPlus2D<T>::BlockPlus2D(BlockF2D<T>& f, BlockF2D<T>& g)
  : BlockCalc2D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockPlus2D<T>::operator()(T output[], const int input[])
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
BlockMinus2D<T>::BlockMinus2D(BlockF2D<T>& f, BlockF2D<T>& g)
  : BlockCalc2D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockMinus2D<T>::operator()(T output[], const int input[])
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
BlockMultiplication2D<T>::BlockMultiplication2D(BlockF2D<T>& f, BlockF2D<T>& g)
  : BlockCalc2D<T>(f,g)
{
  this->getName() = f.getName() + "*" + g.getName();
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockMultiplication2D<T>::operator()(T output[], const int input[])
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
BlockDivision2D<T>::BlockDivision2D (BlockF2D<T>& f, BlockF2D<T>& g)
  : BlockCalc2D<T>(f,g)
{
  this->getName() = f.getName() + "/" + g.getName();
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename T>
bool BlockDivision2D<T>::operator()(T output[], const int input[])
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
BlockF2D<T>& BlockF2D<T>::operator+(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockPlus2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF2D<T>& BlockF2D<T>::operator-(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockMinus2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF2D<T>& BlockF2D<T>::operator*(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockMultiplication2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
BlockF2D<T>& BlockF2D<T>::operator/(BlockF2D<T>& rhs)
{
  auto tmp = std::make_shared< BlockDivision2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}




} // end namespace olb

#endif
