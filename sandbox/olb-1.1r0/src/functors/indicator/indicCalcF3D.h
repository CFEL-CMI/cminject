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

#ifndef INDIC_CALC_F_3D_H
#define INDIC_CALC_F_3D_H


#include "indicatorBaseF3D.h"


namespace olb {


/*
 *  arithmetic helper classes for IndicatorF3D, smoothIndicator3D
 *  UNION         +
 *  WITHOUT       -
 *  INTERSECTION  *
*/

//////////////////////////////// IndicCalc3D ////////////////////////////////
/// arithmetic helper class for Indicator 3d functors
template <typename S>
class IndicCalc3D : public IndicatorF3D<S> {
protected:
  IndicatorF3D<S>& _f;
  IndicatorF3D<S>& _g;
public:
  IndicCalc3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
};

/// addition functor acts as union
template <typename S>
class IndicPlus3D : public IndicCalc3D<S> {
public:
  IndicPlus3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
  bool operator() (bool output[], const S input[]);
};

/// subtraction functor acts as without
template <typename S>
class IndicMinus3D : public IndicCalc3D<S> {
public:
  IndicMinus3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
  bool operator() (bool output[], const S input[]);
};

/// multiplication functor acts as intersection
template <typename S>
class IndicMultiplication3D : public IndicCalc3D<S> {
public:
  IndicMultiplication3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
  bool operator() (bool output[], const S input[]);
};


////////////////////////////////// Discrete Indicator ////////////////////////////////
///// arithmetic helper class for discrete indicator 3d functors
//class DiscIndicCalc3D : public SuperIndicatorF3D {
//protected:
//  SuperIndicatorF3D<T>& _f;
//  SuperIndicatorF3D<T>& _g;
//public:
//  DiscIndicCalc3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g);
//};
//
///// addition functor acts as union
//class DiscIndicPlus3D : public DiscIndicCalc3D {
//public:
//  DiscIndicPlus3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g);
//  bool operator() (bool output[], const int input[]);
//};
//
///// subtraction functor acts as without
//class DiscIndicMinus3D : public DiscIndicCalc3D {
//public:
//  DiscIndicMinus3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g);
//  bool operator() (bool output[], const int input[]);
//};
//
///// multiplication functor acts as intersection
//class DiscIndicMultiplication3D : public DiscIndicCalc3D {
//public:
//  DiscIndicMultiplication3D(SuperIndicatorF3D<T>& f, SuperIndicatorF3D<T>& g);
//  bool operator() (bool output[], const int input[]);
//};



//////////////////////////////// IndicSmoothCalc3D ////////////////////////////////
/// arithmetic helper class for Indicator 3d functors
template <typename T, typename S>
class SmoothIndicCalc3D : public SmoothIndicatorF3D<T,S> {
protected:
  SmoothIndicatorF3D<T,S>& _f;
  SmoothIndicatorF3D<T,S>& _g;
public:
  SmoothIndicCalc3D(SmoothIndicatorF3D<T,S>& f, SmoothIndicatorF3D<T,S>& g);
};

/// addition functor acts as union
template <typename T, typename S>
class SmoothIndicPlus3D : public SmoothIndicCalc3D<T,S> {
public:
  SmoothIndicPlus3D(SmoothIndicatorF3D<T,S>& f, SmoothIndicatorF3D<T,S>& g);
  bool operator() (T output[], const S input[]);
};




} // end namespace olb

#endif
