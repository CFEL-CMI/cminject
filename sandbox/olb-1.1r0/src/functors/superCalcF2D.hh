/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#ifndef SUPER_CALC_F_2D_HH
#define SUPER_CALC_F_2D_HH


#include "functors/superCalcF2D.h"



namespace olb {


template <typename T, typename W>
SuperCalc2D<T,W>::SuperCalc2D(SuperF2D<T,W>& f, SuperF2D<T,W>& g)
  : SuperF2D<T,W>( f.getSuperStructure(), f.getTargetDim() ), _f(f), _g(g)
{
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

// addition
template <typename T, typename W>
SuperPlus2D<T,W>::SuperPlus2D(SuperF2D<T,W>& f, SuperF2D<T,W>& g) : SuperCalc2D<T,W>(f,g)
{
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
}

template <typename T, typename W>
bool SuperPlus2D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]+=tmp[i];
  }
  return true;
}


// subtraction
template <typename T, typename W>
SuperMinus2D<T,W>::SuperMinus2D(SuperF2D<T,W>& f, SuperF2D<T,W>& g) : SuperCalc2D<T,W>(f,g)
{
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T, typename W>
bool SuperMinus2D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]-=tmp[i];
  }
  return true;
}


// multiplication
template <typename T, typename W>
SuperMultiplication2D<T,W>::SuperMultiplication2D(SuperF2D<T,W>& f, SuperF2D<T,W>& g)
  : SuperCalc2D<T,W>(f,g)
{
  this->getName() = f.getName() + "*" + g.getName();
}

template <typename T, typename W>
bool SuperMultiplication2D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]*=tmp[i];
  }
  return true;
}


// division
template <typename T, typename W>
SuperDivision2D<T,W>::SuperDivision2D(SuperF2D<T,W>& f, SuperF2D<T,W>& g)
  : SuperCalc2D<T,W>(f,g)
{
  this->getName() = f.getName() + "/" + g.getName();
}


template <typename T, typename W>
bool SuperDivision2D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]/=tmp[i];
  }
  return true;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator+(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperPlus2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator-(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperMinus2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator*(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperMultiplication2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator/(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperDivision2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}



} // end namespace olb

#endif
