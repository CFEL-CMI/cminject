/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Mathias J. Krause, Cyril Masquelier, Benjamin Förster, Albert Mink
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

#include "indicCalcF3D.h"
#include "indicCalcF3D.hh"


namespace olb {

// arithmetic helper class for indicator 3d functors
template class IndicCalc3D<double>;
template class IndicPlus3D<double>;
template class IndicMinus3D<double>;
template class IndicMultiplication3D<double>;

// arithmetic helper class for Indical 3d functors Smooth
template class SmoothIndicCalc3D<double,double>;
template class SmoothIndicPlus3D<double,double>;

////////////////////// Discrete Indicator //////////////////
//DiscIndicCalc3D::DiscIndicCalc3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g)
//  : _f(f), _g(g)
//{
//  std::swap(f._ptrCalcC, this->_ptrCalcC);
//}
//
//
//DiscIndicPlus3D::DiscIndicPlus3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g)
//  : DiscIndicCalc3D(f, g)
//{}
//
//// returns 1 if( f==1 || g==1 ) UNION
//bool DiscIndicPlus3D::operator() (bool output[], const int input[])
//{
//  this->_f(output, input);
//  bool tmp;
//  this->_g(&tmp, input);
//  output[0] |= tmp;
//  return true;
//}
//
//
//DiscIndicMinus3D::DiscIndicMinus3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g)
//  : DiscIndicCalc3D(f, g)
//{}
//
//// returns 1 if( f==1 && g==0 ) WITHOUT
//bool DiscIndicMinus3D::operator()(bool output[], const int input[])
//{
//  this->_f(output, input);
//  bool tmp;
//  this->_g(&tmp, input);
//  output[0] &= !tmp;
//  return true;
//}
//
//
//// returns 1 if( f==1 && g==1 ) INTERSECTION
//DiscIndicMultiplication3D::DiscIndicMultiplication3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g)
//  : DiscIndicCalc3D(f, g)
//{}
//
//
//bool DiscIndicMultiplication3D::operator() (bool output[], const int input[])
//{
//  this->_f(output, input);
//  bool tmp;
//  this->_g(&tmp, input);
//  output[0] &= tmp;
//  return true;
//}
//
//
//SuperIndicatorF3D<T>& SuperIndicatorF3D::operator+(SuperIndicatorF3D<T>& rhs)
//{
//  auto tmp = std::make_shared< DiscIndicPlus3D >(*this, rhs);
//  this->_ptrCalcC = tmp;
//  return *tmp;
//}
//
//SuperIndicatorF3D<T>& SuperIndicatorF3D::operator-(SuperIndicatorF3D<T>& rhs)
//{
//  auto tmp = std::make_shared< DiscIndicMinus3D >(*this, rhs);
//  this->_ptrCalcC = tmp;
//  return *tmp;
//}
//
//SuperIndicatorF3D<T>& SuperIndicatorF3D::operator*(SuperIndicatorF3D<T>& rhs)
//{
//  auto tmp = std::make_shared< DiscIndicMultiplication3D >(*this, rhs);
//  this->_ptrCalcC = tmp;
//  return *tmp;
//}


}

