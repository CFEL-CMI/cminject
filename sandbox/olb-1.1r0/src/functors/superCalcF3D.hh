/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
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

#ifndef SUPER_CALC_F_3D_HH
#define SUPER_CALC_F_3D_HH


#include "functors/superCalcF3D.h"


namespace olb {


template <typename T, typename W>
SuperCalc3D<T,W>::SuperCalc3D(SuperF3D<T,W>& f, SuperF3D<T,W>& g)
  : SuperF3D<T,W>( f.getSuperStructure(), f.getTargetDim() ), _f(f), _g(g)
{
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

// Addition
template <typename T, typename W>
SuperPlus3D<T,W>::SuperPlus3D(SuperF3D<T,W>& f, SuperF3D<T,W>& g) : SuperCalc3D<T,W>(f,g)
{
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
}

template <typename T, typename W>
bool SuperPlus3D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] += tmp[i];
  }
  return true;
}

// Subtraction
template <typename T, typename W>
SuperMinus3D<T,W>::SuperMinus3D(SuperF3D<T,W>& f, SuperF3D<T,W>& g) : SuperCalc3D<T,W>(f,g)
{
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T, typename W>
bool SuperMinus3D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] -= tmp[i];
  }
  return true;
}


// Subtraction (bool sepcialzation)
template <typename T>
SuperMinus3D<T,bool>::SuperMinus3D(SuperF3D<T,bool>& f, SuperF3D<T,bool>& g) : SuperCalc3D<T,bool>(f,g)
{
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
}
// Subtraction specialization for bool -> "Without" operator
template <typename T>
bool SuperMinus3D<T,bool>::operator()(bool output[], const int input[])
{
  this->_f(output,input);
  bool tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] &= !tmp[i];
  }
  return true;
}


// Multiplication
template <typename T, typename W>
SuperMultiplication3D<T,W>::SuperMultiplication3D(SuperF3D<T,W>& f, SuperF3D<T,W>& g)
  : SuperCalc3D<T,W>(f,g)
{
  this->getName() = f.getName() + "*" + g.getName();
}

template <typename T, typename W>
bool SuperMultiplication3D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] *= tmp[i];
  }
  return true;
}


// Division
template <typename T, typename W>
SuperDivision3D<T,W>::SuperDivision3D(SuperF3D<T,W>& f, SuperF3D<T,W>& g)
  : SuperCalc3D<T,W>(f,g)
{
  this->getName() = f.getName() + "/" + g.getName();
}

template <typename T, typename W>
bool SuperDivision3D<T,W>::operator()(W output[], const int input[])
{
  this->_f(output,input);
  W tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] /= tmp[i];
  }
  return true;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator+(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperPlus3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator-(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperMinus3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator*(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperMultiplication3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF3D<T,W>& SuperF3D<T,W>::operator/(SuperF3D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperDivision3D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}



} // end namespace olb

#endif
