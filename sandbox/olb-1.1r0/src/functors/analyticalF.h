/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_F_H
#define ANALYTICAL_F_H


#include<vector>
#include "functors/analyticalBaseF.h"
#include "functors/indicator/indicatorF2D.h"
#include "functors/indicator/indicatorF3D.h"

/**
 *  The functor dimensions are given by F: S^m -> T^n  (S=source, T=target)
 *  and are implemented via GenericF(n,m).
 *  Don't get confused by the flipped order of source and target.
 */

namespace olb {

template<typename T, typename S> class SmoothIndicatorSphere3D;

////////////////////////////////////////////////////////////////////////////////
////////implementation of several 1d,2d,3d functors (analyticalFXD)/////////////
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////1D////////////////////////////////////////////

/// AnalyticalConst1D: 1D -> XD, where XD is defined by value.size()
template <typename T, typename S>
class AnalyticalConst1D : public AnalyticalF1D<T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _c;
public:
  AnalyticalConst1D(T value);
  AnalyticalConst1D(std::vector<T>& value);
  bool operator() (T output[], const S x[]);
};

/// AnalyticalLinear1D: 1D -> 1D troughout given points (x0,v0) and (x1,v1)
//  Punktsteigungsform
template <typename T, typename S>
class AnalyticalLinear1D : public AnalyticalF1D<T,S> {
private:
  T _a;
  T _b;
public:
  AnalyticalLinear1D(T a, T b);
  AnalyticalLinear1D(S x0, T v0, S x1, T v1);
  bool operator() (T output[], const S x[]); ///< returns line _a*x + _b
};

/// AnalyticalRandom1D: 1D -> 1D with random image in (0,1)
template <typename T, typename S>
class AnalyticalRandom1D : public AnalyticalF1D<T,S> {
public:
  AnalyticalRandom1D();
  bool operator() (T output[], const S x[]);
};


/// represents an inverse parabola profile like it is used in Poiseuille inflow
/// note: output depends only on first parameter, maps 1D,2D,3D->1D
template <typename T, typename S>
class AnalyticalSquare1D : public AnalyticalF1D<T,S> {
private:
  S _cp;
  S _r;
  T _maxi;
public:
  AnalyticalSquare1D(S cp, S r, T maxi);
  bool operator() (T output[], const S x[]);
};


/// SinusStartScale: 1D -> 1D a start curve based on sinus for a continuous transition at 0 and 1
template <typename T, typename S>
class SinusStartScale : public AnalyticalF1D<T,S> {
protected:
  S _numTimeSteps;
  T _maxValue;
  T _pi;
public:
  SinusStartScale(int numTimeSteps=1, T maxValue=1);
  bool operator() (T output[], const S x[]);
};


/// PolynomialStartScale: 1D -> 1D a start curve based on a polynomial fifth order for a continuous transition at 0 and 1: maxValue*(6*y^5-15*y^4+10*y^3)
template <typename T, typename S>
class PolynomialStartScale : public AnalyticalF1D<T,S> {
protected:
  S _numTimeSteps;
  T _maxValue;
public:
  PolynomialStartScale(S numTimeSteps=S(1), T maxValue=T(1));
  bool operator() (T output[], const S x[]);
};

/// Derivative of a given 1D functor computed with a finite difference
template <typename T>
class AnalyticalDiffFD1D : public AnalyticalF1D<T,T> {
protected:
  AnalyticalF1D<T,T>& _f;
  T _eps;
public:
  AnalyticalDiffFD1D(AnalyticalF1D<T,T>& f, T eps = 1.e-10);
  bool operator() (T output[], const T input[]);
};

//////////////////////////////////2D////////////////////////////////////////////

template <typename T, typename S>
class AnalyticalComposed2D final : public AnalyticalF2D<T,S> {
private:
  AnalyticalF2D<T,S>& _f0;
  AnalyticalF2D<T,S>& _f1;
public:
  AnalyticalComposed2D(AnalyticalF2D<T,S>& f0, AnalyticalF2D<T,S>& f1);
  bool operator() (T output[], const S x[]);
};

/// AnalyticalConst2D: 2D -> XD, where XD is defined by value.size()
template <typename T, typename S>
class AnalyticalConst2D final : public AnalyticalF2D<T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _c;
public:
  AnalyticalConst2D(T value);
  AnalyticalConst2D(T value0, T value1);
  AnalyticalConst2D(std::vector<T>& value);
  bool operator() (T output[], const S x[]);
};

/// AnalyticalLinear2D: 2D -> 1D troughout given points (x0,y0,v0), (x1,y1,v1), (x2,y2,v2)
template <typename T, typename S>
class AnalyticalLinear2D final : public AnalyticalF2D<T,S> {
protected:
  T _a;
  T _b;
  T _c;
public:
  AnalyticalLinear2D(T a, T b, T c);
  AnalyticalLinear2D(S x0, S y0, T v0, S x1, S y1, T v1, S x2, S y2, T v2);
  bool operator() (T output[], const S x[]);
};

/// AnalyticalRandom2D: 2D -> 1D with random image in (0,1)
template <typename T, typename S>
class AnalyticalRandom2D final : public AnalyticalF2D<T,S> {
public:
  AnalyticalRandom2D();
  bool operator() (T output[], const S x[]);
};

/// AnalyticalRandom2D: 2D -> 1D with maxValue in the center decreasing linearly with the distrance to the center to zero at the radius and zero outside
template <typename T, typename S>
class AnalyticalParticleAdsorptionLinear2D final : public AnalyticalF2D<T,S> {
protected:
  T _center[2];
  T _radius;
  T _maxValue;
public:
  AnalyticalParticleAdsorptionLinear2D(T center[], T radius, T maxValue);
  bool operator() (T output[], const S x[]);
};

template <typename T, typename S>
class ParticleU2D : public AnalyticalF2D<T,S> {
protected:
  ParticleIndicatorF2D<T,T>& _indicator;
  std::vector<T> _u;
  T _omega;
public:
  ParticleU2D(ParticleIndicatorF2D<T,T>& indicator, std::vector<T>& u, T& omega);
  bool operator()(T output[], const S input[]);
};


//////////////////////////////////3D////////////////////////////////////////////
template <typename T, typename S>
class AnalyticalComposed3D final : public AnalyticalF3D<T,S> {
private:
  AnalyticalF3D<T,S>& _f0;
  AnalyticalF3D<T,S>& _f1;
  AnalyticalF3D<T,S>& _f2;
public:
  AnalyticalComposed3D(AnalyticalF3D<T,S>& f0, AnalyticalF3D<T,S>& f1, AnalyticalF3D<T,S>& f2);
  bool operator() (T output[], const S x[]);
};

/// AnalyticalConst3D: 3D -> XD, where XD is defined by value.size()
template <typename T, typename S>
class AnalyticalConst3D final : public AnalyticalF3D<T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _c;
public:
  AnalyticalConst3D(T value);
  AnalyticalConst3D(T value0, T value1);
  AnalyticalConst3D(T value0, T value1, T value2);
  AnalyticalConst3D(std::vector<T>& value);
  bool operator() (T output[], const S x[]);
};

/// AnalyticalLinear3D: 3D -> 1D troughout given points (x0,y0,z0,v0), (x1,y1,z1,v1), (x2,y2,z2,v2), (x3,y3,z3,v3)
template <typename T, typename S>
class AnalyticalLinear3D final : public AnalyticalF3D<T,S> {
protected:
  T _a;
  T _b;
  T _c;
  T _d;
public:
  AnalyticalLinear3D(T a, T b, T c, T d);
  AnalyticalLinear3D(S x0, S y0, S z0, T v0, S x1, S y1, S z1, T v1, S x2, S y2,
                     S z2, T v2, S x3, S y3, S z3, T v3);
  bool operator() (T output[], const S x[]);
};

/// AnalyticalRandom3D: 3D -> 1D with random image in (0,1)
template <typename T, typename S>
class AnalyticalRandom3D final : public AnalyticalF3D<T,S> {
public:
  AnalyticalRandom3D();
  bool operator() (T output[], const S x[]);
};

/// AnalyticalScaled3D: 3D -> Image(AnalyticalF) scales AnalyticalF by _scale
template <typename T, typename S>
class AnalyticalScaled3D final : public AnalyticalF3D<T,S> {
private:
  AnalyticalF3D<T,S>& _f;
  T _scale;
public:
  AnalyticalScaled3D(AnalyticalF3D<T,S>& f, T scale);
  bool operator() (T output[], const S x[]);
};

/// TODO Comment
template <typename T, typename S>
class ParticleU3D : public AnalyticalF3D<T,S> {
protected:
  ParticleIndicatorF3D<T,T>& _indicator;
  std::vector<T> _u;
  std::vector<T> _omega;
public:
  ParticleU3D(ParticleIndicatorF3D<T,T>& indicator, std::vector<T>& u, std::vector<T>& omega);
  bool operator()(T output[], const S input[]);
};















// some ideas...

/// represents bilinear functions from 1D -> 2D
template <typename T, typename S>
class BilinearAnalyticalF : public AnalyticalF1D<T,S> { };


/*
/// represents arbitrary polynomial functions
template <typename T, typename S>
class PolynomialAnalyticalF : public AnalyticalF1D<T,S> {
private:
  std::vector<T> c;
public:
  PolynomialAnalyticalF(std::vector<int> coeffs) : c(coeffs), AnalyticalF1D<T,S>(2,1) { }
  T operator() (S x) {
    T sum = 0;
    for (int i=0; i<size(coeffs); i++) {
      sum += coeffs(i)*pow(x,i);                    // sum+= c(i)*x^i
    }
    return sum;
  }
};
*/



} // end namespace olb

#endif
