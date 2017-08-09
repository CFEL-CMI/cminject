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

#ifndef ANALYTICAL_F_HH
#define ANALYTICAL_F_HH

#include<vector>
#include<cmath>     // for lpnorm
#include<stdlib.h>  // for random
#include <iostream>
#include<time.h>

#include "functors/analyticalF.h"
#include "core/singleton.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {


//////////////////////////////////1D////////////////////////////////////////////
template <typename T, typename S>
AnalyticalConst1D<T,S>::AnalyticalConst1D(std::vector<T>& value)
  : AnalyticalF1D<T,S>(value.size()), _c(value)
{
  this->getName() = "const";
}

template <typename T, typename S>
AnalyticalConst1D<T,S>::AnalyticalConst1D(T value) : AnalyticalF1D<T,S>(1)
{
  _c.push_back(value);
  this->getName() = "const";
}

template <typename T, typename S>
bool AnalyticalConst1D<T,S>::operator()(T output[], const S x[])
{
  for ( unsigned i = 0; i < _c.size(); ++i) {
    output[i] = _c[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticalLinear1D<T,S>::AnalyticalLinear1D(T a, T b)
  : AnalyticalF1D<T,S>(1), _a(a), _b(b)
{
  this->getName() = "linear";
}

template <typename T, typename S>
AnalyticalLinear1D<T,S>::AnalyticalLinear1D(S x0, T v0, S x1, T v1)
  : AnalyticalF1D<T,S>(1)
{
  if ( util::nearZero(x1-x0) ) {
    std::cout << "Error: x1-x2=0" << std::endl;
  } else {
    _a = ( v1-v0 ) / ( x1-x0 );
    _b = v0 - _a*x0;
  }
  this->getName() = "linear";
}

template <typename T, typename S>
bool AnalyticalLinear1D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_a*x[0] + _b;
  return true;
}


template <typename T, typename S>
AnalyticalRandom1D<T,S>::AnalyticalRandom1D() : AnalyticalF1D<T,S>(1)
{
#ifdef PARALLEL_MODE_MPI
  int  nameLen, numProcs, myID;
  char processorName[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_rank(MPI_COMM_WORLD,&myID);
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Get_processor_name(processorName,&nameLen);
  srand(time(0)+myID*numProcs + nameLen);
#else
  srand(time(0));
#endif
  this->getName() = "random";
}

template <typename T, typename S>
bool AnalyticalRandom1D<T,S>::operator()(T output[], const S x[])
{
  output[0]=(rand()%RAND_MAX)/(T)RAND_MAX;
  return true;
}


template <typename T, typename S>
AnalyticalSquare1D<T,S>::AnalyticalSquare1D(S cp, S r, T maxi)
  : AnalyticalF1D<T,S>(1), _cp(cp), _r(r), _maxi(maxi)
{
  this->getName() = "square";
}

template <typename T, typename S>
bool AnalyticalSquare1D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_maxi*(1-(x[0]-_cp) * (x[0]-_cp) / (T)_r / (T)_r);
  return true;
}



/////////////////////////////////someOtherFunctors//////////////////////////////
template <typename T, typename S>
PolynomialStartScale<T,S>::PolynomialStartScale(S numTimeSteps, T maxValue)
  : AnalyticalF1D<T,S>(1), _numTimeSteps(numTimeSteps), _maxValue(maxValue)
{
  this->getName() = "polyStartScale";
}

template <typename T, typename S>
bool PolynomialStartScale<T,S>::operator()(T output[], const S x[])
{
  output[0]=(T) x[0] / (T) _numTimeSteps;
  output[0]=_maxValue * output[0]*output[0]*output[0] * ( 10.0 + output[0] * (6.0*output[0] - 15.0) );
  return true;
}


template <typename T, typename S>
SinusStartScale<T,S>::SinusStartScale(int numTimeSteps, T maxValue)
  : AnalyticalF1D<T,S>(1), _numTimeSteps(numTimeSteps), _maxValue(maxValue),
    _pi(4.0*atan(1.0))
{
  this->getName() = "sinusStartScale";
}

template <typename T, typename S>
bool SinusStartScale<T,S>::operator()(T output[], const S x[])
{
  output[0]=(_maxValue * (sin(-_pi / 2.0 + (T)x[0] / (T)_numTimeSteps * _pi) + 1.0)) / 2.0;
  return true;
}

template <typename T>
AnalyticalDiffFD1D<T>::AnalyticalDiffFD1D(AnalyticalF1D<T,T>& f, T eps) : AnalyticalF1D<T,T>(1), _f(f), _eps(eps)
{
}

template <typename T>
bool AnalyticalDiffFD1D<T>::operator() (T output[], const T input[])
{
  _f(output,input);
  T x = output[0];
  T input2[1];
  input2[0] = input[0] + _eps;
  _f(output,input2);
  output[0] -= x;
  output[0] /= _eps;
  return true;
}

///////////////////////////////////////2D///////////////////////////////////////
template <typename T, typename S>
AnalyticalComposed2D<T,S>::AnalyticalComposed2D( AnalyticalF2D<T,S>& f0,
    AnalyticalF2D<T,S>& f1)
  : AnalyticalF2D<T,S>(2), _f0(f0), _f1(f1)
{
  this->getName() = "composed";
}

template <typename T, typename S>
bool AnalyticalComposed2D<T,S>::operator()(T output[], const S x[])
{
  T tmp1[_f0.getTargetDim()], tmp2[_f1.getTargetDim()];
  _f0(tmp1,x);
  _f1(tmp2,x);
  output[0]=tmp1[0];
  output[1]=tmp2[0];
  return true;
}


template <typename T, typename S>
AnalyticalConst2D<T,S>::AnalyticalConst2D(std::vector<T>& value)
  : AnalyticalF2D<T,S>(value.size()), _c(value)
{
  this->getName() = "const";
}

template <typename T, typename S>
AnalyticalConst2D<T,S>::AnalyticalConst2D(T value) : AnalyticalF2D<T,S>(1)
{
  _c.push_back(value);
  this->getName() = "const";
}

template <typename T, typename S>
AnalyticalConst2D<T,S>::AnalyticalConst2D(T value0, T value1) : AnalyticalF2D<T,S>(2)
{
  _c.push_back(value0);
  _c.push_back(value1);
  this->getName() = "const";
}

template <typename T, typename S>
bool AnalyticalConst2D<T,S>::operator()(T output[], const S x[])
{
  for (unsigned i = 0; i < _c.size(); ++i) {
    output[i] = _c[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticalLinear2D<T,S>::AnalyticalLinear2D(T a, T b, T c)
  : AnalyticalF2D<T,S>(1), _a(a), _b(b), _c(c)
{
  this->getName() = "linear";
}

template <typename T, typename S>
AnalyticalLinear2D<T,S>::AnalyticalLinear2D(S x0, S y0, T v0, S x1, S y1,
    T v1, S x2, S y2, T v2)
  : AnalyticalF2D<T,S>(1)
{
  this->getName() = "linear";
  T n2= (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
  if ( util::nearZero(n2) ) {
    std::cout << "Error function" << std::endl;
  } else {
    T n0 = (y1-y0)*(v2-v0) - (v1-v0)*(y2-y0);
    T n1 = (v1-v0)*(x2-x0) - (x1-x0)*(v2-v0);
    _a = -n0 / n2;
    _b = -n1 / n2;
    _c = (x0*n0 + y0*n1 + v0*n2) / n2;
  }
}

template <typename T, typename S>
bool AnalyticalLinear2D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_a*x[0] + _b*x[1] + _c;
  return true;
}

template <typename T, typename S>
AnalyticalParticleAdsorptionLinear2D<T,S>::AnalyticalParticleAdsorptionLinear2D(T center[], T radius, T maxValue) : AnalyticalF2D<T,S>(2), _radius(radius), _maxValue(maxValue)
{
  _center[0] = center[0];
  _center[1] = center[1];
  this->getName() = "particleAdsorptionLinear2D";
}

template <typename T, typename S>
bool AnalyticalParticleAdsorptionLinear2D<T,S>::operator()(T output[], const S input[])
{
  T dist = sqrt((input[0]-_center[0])*(input[0]-_center[0]) + (input[1]-_center[1])*(input[1]-_center[1]));

  if (dist > _radius) {
    output[0] = T();
    output[1] = T();
    return true;
  } else {
    output[0] = _maxValue*(T(1) - dist/_radius)*(_center[0]-input[0])/_radius;
    output[1] = _maxValue*(T(1) - dist/_radius)*(_center[1]-input[1])/_radius;
    return true;
  }
}



template <typename T, typename S>
AnalyticalRandom2D<T,S>::AnalyticalRandom2D() : AnalyticalF2D<T,S>(1)
{
#ifdef PARALLEL_MODE_MPI
  int  nameLen, numProcs, myID;
  char processorName[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_rank(MPI_COMM_WORLD,&myID);
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Get_processor_name(processorName,&nameLen);
  srand(time(0)+myID*numProcs + nameLen);
#else
  srand(time(0));
#endif
  this->getName() = "random";
}

template <typename T, typename S>
bool AnalyticalRandom2D<T,S>::operator()(T output[], const S x[])
{
  output[0]=(rand()%RAND_MAX)/(T)RAND_MAX;
  return true;
}

template <typename T, typename S>
ParticleU2D<T,S>::ParticleU2D(ParticleIndicatorF2D<T,T>& indicator, std::vector<T>& u, T& omega)
  :AnalyticalF2D<T,S>(2), _indicator(indicator), _u(u), _omega(omega)
{
  this->getName() = "ParticleU";
}

template <typename T, typename S>
bool ParticleU2D<T,S>::operator()(T output[], const S input[])
{
  //T inside[1];
  output[0] = T();
  output[1] = T();

  //_indicator(inside, input);
  //if (inside[0] != 0) {

  //two dimensions: u = U + w x r = (Ux, Uy, 0) + (0,0,w) x (x,y,0) = (Ux, Uy, 0) + (-wy, wx, 0)
  output[0] = _u[0] - _omega*(input[1] - _indicator.getPos()[1]);
  output[1] = _u[1] + _omega*(input[0] - _indicator.getPos()[0]);

  return true;
}


///////////////////////////////////////3D///////////////////////////////////////
template <typename T, typename S>
AnalyticalComposed3D<T,S>::AnalyticalComposed3D(AnalyticalF3D<T,S>& f0,
    AnalyticalF3D<T,S>& f1, AnalyticalF3D<T,S>& f2)
  : AnalyticalF3D<T,S>(3), _f0(f0), _f1(f1), _f2(f2)
{
  this->getName() = "composed";
}

template <typename T, typename S>
bool AnalyticalComposed3D<T,S>::operator()(T output[], const S x[])
{
  T outputTmp[_f0.getTargetDim()];
  _f0(outputTmp,x);
  output[0]=outputTmp[0];
  _f1(outputTmp,x);
  output[1]=outputTmp[0];
  _f2(outputTmp,x);
  output[2]=outputTmp[0];
  return true;
}


template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(std::vector<T>& value)
  : AnalyticalF3D<T,S>(value.size()), _c(value)
{
  this->getName() = "const";
}

template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(T value) : AnalyticalF3D<T,S>(1)
{
  _c.push_back(value);
  this->getName() = "const";
}

template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(T value0, T value1) : AnalyticalF3D<T,S>(2)
{
  _c.push_back(value0);
  _c.push_back(value1);
  this->getName() = "const";
}

template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(T value0, T value1, T value2)
  : AnalyticalF3D<T,S>(3)
{
  _c.push_back(value0);
  _c.push_back(value1);
  _c.push_back(value2);
  this->getName() = "const";
}

template <typename T, typename S>
bool AnalyticalConst3D<T,S>::operator()(T output[], const S x[])
{
  for (unsigned int i = 0; i < _c.size(); ++i) {
    output[i] = _c[i];
  }
  return true;
}


template <typename T, typename S>
AnalyticalLinear3D<T,S>::AnalyticalLinear3D(T a, T b, T c, T d)
  : AnalyticalF3D<T,S>(1), _a(a), _b(b), _c(c), _d(d)
{
  this->getName() = "linear";
}

template <typename T, typename S>
AnalyticalLinear3D<T,S>::AnalyticalLinear3D(S x0, S y0, S z0, T v0, S x1,
    S y1, S z1, T v1, S x2, S y2, S z2, T v2, S x3, S y3, S z3, T v3)
  : AnalyticalF3D<T,S>(1)
{
  this->getName() = "linear";
  T n = ( (y3-y0)*(x1-x0)-(x3-x0)*(y1-y0) ) * ( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) )
        +( (z3-z0)*(x1-x0)-(x3-x0)*(z1-z0) ) * ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
  if ( util::nearZero(n) ) {
    std::cout << "Error function" << std::endl;
  } else {
    T w = ( (y1-y0)*(x3-x0)-(x1-x0)*(y3-y0) ) * ( (v2-v0)-(x2-x0)*(v1-v0) / (x1-x0) )
          /( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) ) + (v3-v0) - (x3-x0)*(v1-v0) / (x1-x0);
    T zx = (y1-y0)*( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) )
           -(z1-z0)*( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T rx = (v1-v0)/(x1-x0) - (y1-y0)*(v2-v0) / ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) )
           +(y1-y0)*(x2-x0)*(v1-v0) / ( (y2-y0)*(x1-x0)*(x1-x0)-(x2-x0)*(y1-y0)*(x1-x0) );
    T zy = (x1-x0)*( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) );
    T ry = ( (x1-x0)*(v2-v0)-(x2-x0)*(v1-v0) ) / ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T zz = (x1-x0)*( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T h = w/n;
    _a = rx + zx*h;
    _b = ry + zy*h;
    _c = zz*h;
    _d = v0 - x0*_a - y0*_b - z0*_c;
  }
}

template <typename T, typename S>
bool AnalyticalLinear3D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_a*x[0] + _b*x[1] + _c*x[2] + _d;
  return true;
}


template <typename T, typename S>
AnalyticalRandom3D<T,S>::AnalyticalRandom3D() : AnalyticalF3D<T,S>(1)
{
#ifdef PARALLEL_MODE_MPI
  int  nameLen, numProcs, myID;
  char processorName[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_rank(MPI_COMM_WORLD,&myID);
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Get_processor_name(processorName,&nameLen);
  srand(time(0)+myID*numProcs + nameLen);
#else
  srand(time(0));
#endif
  this->getName() = "random";
}

template <typename T, typename S>
bool AnalyticalRandom3D<T,S>::operator()(T output[], const S x[])
{
  output[0]=(rand()%RAND_MAX)/(T)RAND_MAX;
  return true;
}


template <typename T, typename S>
AnalyticalScaled3D<T,S>::AnalyticalScaled3D(AnalyticalF3D<T,S>& f, T scale)
  : AnalyticalF3D<T,S>(f.getTargetDim()), _f(f), _scale(scale)
{
  this->getName() = "scaled";
}

template <typename T, typename S>
bool AnalyticalScaled3D<T,S>::operator()(T output[], const S x[])
{
  T outputTmp[_f.getTargetDim()];
  _f(outputTmp,x);
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = _scale*outputTmp[iDim];
  }
  return true;
}


template <typename T, typename S>
ParticleU3D<T,S>::ParticleU3D(ParticleIndicatorF3D<T,T>& indicator,
                              std::vector<T>& u, std::vector<T>& omega)
  : AnalyticalF3D<T,S>(1), _indicator(indicator), _u(u), _omega(omega)
{
  this->getName() = "ParticleU";
}


template <typename T, typename S>
bool ParticleU3D<T,S>::operator()(T output[], const S input[])
{
  if (_indicator(output, input) != 0) {
    output[0] = _u[0] + ( _omega[1]*(input[2] - _indicator.getPos()[2]) - _omega[2]*(input[1] - _indicator.getPos()[1]) );
    output[1] = _u[1] + ( _omega[2]*(input[0] - _indicator.getPos()[0]) - _omega[0]*(input[2] - _indicator.getPos()[2]) );
    output[2] = _u[2] + ( _omega[0]*(input[1] - _indicator.getPos()[1]) - _omega[1]*(input[0] - _indicator.getPos()[0]) );
  }
  return true;
}




} // end namespace olb

#endif
