/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Cyril Masquelier, Jan Marquardt, Mathias J. Krause
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

#ifndef INDICATOR_F_2D_HH
#define INDICATOR_F_2D_HH


#include<cmath>
#include "indicatorF2D.h"
#include "utilities/vectorHelpers.h"
#include "functors/interpolationF3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define M_PI2 1.57079632679489661923

namespace olb {

// Warning : the cuboid is only defined parallel to the plans x=0 and y=0 !!!
// For other cuboids, please use the parallelepiped version
template <typename S>
IndicatorCuboid2D<S>::IndicatorCuboid2D(Vector<S,2> extend, Vector<S,2> origin, S theta)
  : _theta(theta)
{
  this->_myMin = origin;
  this->_myMax = origin + extend;
  _center = origin + S(.5)*extend;
  _xLength = extend[0];
  _yLength = extend[1];
}

template <typename S>
IndicatorCuboid2D<S>::IndicatorCuboid2D(S xLength, S yLength, Vector<S,2> center, S theta )
  : _center(center), _xLength(xLength), _yLength(yLength), _theta(-theta)
{
  this->_myMin = {_center[0] - _xLength/2., _center[1] - _yLength/2.};
  this->_myMax = {_center[0] + _xLength/2., _center[1] + _yLength/2.};
}


// returns true if x is inside the cuboid
template <typename S>
bool IndicatorCuboid2D<S>::operator()(bool output[], const S input[])
{
  S x, y;
  if ( !util::nearZero(_theta) ) {
    x = _center[0] + (input[0] - _center[0])*std::cos(_theta) - (input[1] - _center[1])*std::sin(_theta);
    y = _center[1] + (input[0] - _center[0])*std::sin(_theta) + (input[1] - _center[1])*std::cos(_theta);
  } else {
    x = input[0];
    y = input[1];
  }

  output[0] = (  (fabs(_center[0] - x) < _xLength/2. || nearZero(fabs(_center[0] - x) - _xLength/2.))
                 && (fabs(_center[1] - y) < _yLength/2. || nearZero(fabs(_center[1] - y) - _yLength/2.)) );
  return true;
}

// creator function
template <typename S>
IndicatorCuboid2D<S>* createIndicatorCuboid2D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCuboid2D");

  Vector<S,2> origin;
  Vector<S,2> extend(1);

  params.setWarningsOn(false);
  // params[parameter].read(value) sets the value or returns false if the parameter can not be found

  if (!params["originX"].read(origin[0], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: originX. Set default: originX=0" << std::endl;
  }
  if (!params["originY"].read(origin[1], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: originY. Set default: originY=0" << std::endl;
  }


  if (!params["extendX"].read(extend[0], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: extendX. Set default: extendX=1" << std::endl;
  }
  if (!params["extendY"].read(extend[1], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: extendY. Set default: extendY=1" << std::endl;
  }

  params.setWarningsOn(true);

  return new IndicatorCuboid2D<S>(extend, origin);
}


template <typename S>
IndicatorCircle2D<S>::IndicatorCircle2D(Vector<S,2> center, S radius)
  :  _center(center),
     _radius2(radius*radius)
{
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

// returns true if x is inside the circle
template <typename S>
bool IndicatorCircle2D<S>::operator()(bool output[], const S input[])
{
  output[0] = ( std::pow(_center[0] - input[0],2) + std::pow(_center[1] - input[1], 2) <= _radius2 );
  return output[0];
}


template <typename S>
bool IndicatorCircle2D<S>::distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC)
{
  S a = direction[0]*direction[0] + direction[1]*direction[1];

  // returns 0 if point is at the boundary of the sphere
  if ( util::nearZero(a-_radius2) ) {
    distance = S(0);
    return true;
  }
  // norm of direction
  a = sqrt(a);

  S b = 2.*((origin[0] - _center[0])*direction[0] +
            (origin[1] - _center[1])*direction[1])/a;
  S c = -_radius2 + (origin[0] - _center[0])*(origin[0] - _center[0])
        + (origin[1] - _center[1])*(origin[1] - _center[1]);

  // discriminant
  S d = b*b - 4.*c;
  if (d < 0) {
    return false;
  }

  S x1 = (- b + sqrt(d)) *0.5;
  S x2 = (- b - sqrt(d)) *0.5;

  // case if origin is inside the sphere
  if ((x1<0.) || (x2<0.)) {
    if (x1>0.) {
      distance = x1;
      return true;
    }
    if (x2>0.) {
      distance = x2;
      return true;
    }
  }
  // case if origin is ouside the sphere
  else {
    distance = std::min(x1,x2);
    return true;
  }

  return false;
}


template <typename S>
bool IndicatorCircle2D<S>::normal(Vector<S,2>& normal, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC)
{
  S dist;
  if (!(distance(dist, origin, direction, iC)) ) {
    return false;
  }

  Vector<S,2> intresection(origin + dist*direction);

  normal = intresection - _center;

  return true;
}



template <typename S>
IndicatorCircle2D<S>* createIndicatorCircle2D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCircle2D");

  Vector<S,2> center;
  S radius = 1;

  params.setWarningsOn(false);
  // params[parameter].read(value) sets the value or returns false if the parameter can not be found

  if (!params["Application"]["circle"]["centerX"].read(center[0], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: centerX. Set default: centerX=0" << std::endl;
  }
  if (!params["Application"]["circle"]["centerY"].read(center[1], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: centerY. Set default: centerY=0" << std::endl;
  }

  if (!params["Application"]["circle"]["radius"].read(radius, verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: radius. Set default: radius=1" << std::endl;
  }


  return new IndicatorCircle2D<S>(center, radius);
}


template <typename T, typename S>
SmoothIndicatorCuboid2D<T,S>::SmoothIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S mass, S epsilon, S theta)
  : _xLength(xLength),_yLength(yLength)
{
  this->_center = center;
  this->_epsilon = epsilon;
  this->_theta = theta;
  this->_mass = mass;
  this->_mofi = 1./12.*this->_mass*(_xLength*_xLength+_yLength*_yLength);
  this->_myMin = {this->_center[0] - getRadius() - 10*this->_epsilon, this->_center[1] - getRadius() - 10*this->_epsilon};
  this->_myMax = {this->_center[0] + getRadius() + 10*this->_epsilon, this->_center[1] + getRadius() + 10*this->_epsilon};
}

template <typename T, typename S>
bool SmoothIndicatorCuboid2D<T,S>::operator()(T output[], const S r[])
{
  T xDist = r[0] - this->_center[0];
  T yDist = r[1] - this->_center[1];

  T xL2 = _xLength/2.;
  T yL2 = _yLength/2.;

  // counter-clockwise rotation by _theta=-theta around center
  T ct = std::cos(this->_theta);
  T st = std::sin(this->_theta);

  T x = this->_center[0] + xDist*ct - yDist*st;
  T y = this->_center[1] + xDist*st + yDist*ct;

  xDist = fabs(x -this-> _center[0]);
  yDist = fabs(y -this-> _center[1]);

  if ( xDist <= xL2 && yDist <= yL2) {
    output[0] = 1.;
    return true;
  }
  if ( xDist > xL2 + this->_epsilon || yDist > yL2 + this->_epsilon ) {
    output[0] = 0.;
    return false;
  }
  if ( xDist < xL2 && (yDist <= yL2 + this->_epsilon  && yDist > yL2) ) {
    output[0] = T( std::pow(cos(M_PI2*(yDist - yL2)/this->_epsilon), 2));
    return true;
  }
  if ( yDist < yL2 && (xDist <= xL2 + this->_epsilon  && xDist > xL2) ) {
    output[0] = T( std::pow(cos(M_PI2*(xDist - xL2)/this->_epsilon), 2));
    return true;
  }
  if ( (xDist <= xL2 + this->_epsilon && xDist > xL2) && (yDist <= yL2 + this->_epsilon && yDist > yL2) ) {
    output[0] = T( (std::pow(cos(M_PI2*(xDist - xL2)/this->_epsilon), 2) *
                    std::pow(cos(M_PI2*(yDist - yL2)/this->_epsilon), 2)) );
    return true;
  }
  output[0] = 0.;
  return false;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorCuboid2D<T,S>::getMin()
{
  this->_myMin[0] = this->_center[0] - getRadius() - 10*this->_epsilon;
  this->_myMin[1] = this->_center[1] - getRadius() - 10*this->_epsilon;
  return this->_myMin;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorCuboid2D<T,S>::getMax()
{
  this->_myMax[0] = this->_center[0] + getRadius() + 10*this->_epsilon;
  this->_myMax[1] = this->_center[1] + getRadius() + 10*this->_epsilon;
  return this->_myMax;
}

template <typename T, typename S>
S SmoothIndicatorCuboid2D<T,S>::getRadius()
{
  return .5*(std::sqrt(std::pow(_xLength+this->_epsilon, 2)+std::pow(_yLength+this->_epsilon, 2)));
}

template <typename T, typename S>
S SmoothIndicatorCuboid2D<T,S>::getDiam()
{
  return (std::sqrt(std::pow(_xLength, 2)+std::pow(_yLength, 2)));
}


template <typename T, typename S>
SmoothIndicatorCircle2D<T,S>::SmoothIndicatorCircle2D(Vector<S,2> center, S radius, S mass, S epsilon)
{
  this->_radius = radius;
  this->_epsilon = epsilon;
  this->_center = center;
  this->_myMin = this->_center - this->_radius - this->_epsilon *0.5;
  this->_myMax = this->_center + this->_radius + this->_epsilon *0.5;
  this->_mass = mass;
  this->_mofi = 0.5 * this->_mass * pow(this->_radius, 2);
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorCircle2D<T,S>::getMin()
{
  this->_myMin[0] = this->_center[0] - this->_radius - .5*this->_epsilon;
  this->_myMin[1] = this->_center[1] - this->_radius - .5*this->_epsilon;
  return this->_myMin;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorCircle2D<T,S>::getMax()
{
  this->_myMax[0] = this->_center[0] + this->_radius + .5*this->_epsilon;
  this->_myMax[1] = this->_center[1] + this->_radius + .5*this->_epsilon;
  return this->_myMax;
}

// returns true if x is inside the sphere
template <typename T, typename S>
bool SmoothIndicatorCircle2D<T,S>::operator()(T output[], const S input[])
{
  double d;   // distance to the figure
  double distToCenter2 = std::pow(this->_center[0]-input[0], 2) +
                         std::pow(this->_center[1]-input[1], 2);
  if ( distToCenter2 >= std::pow(this->_radius + this->_epsilon *0.5, 2)) {
    output[0] = T(0);
    return true;
  } else if ( distToCenter2 <= std::pow(this->_radius - this->_epsilon *0.5, 2)) {
    output[0] = T(1);
    return true;
  } else {
    // d is between 0 and _epsilon
    d = std::sqrt(distToCenter2) - this->_radius + this->_epsilon *0.5;
    output[0] = T(std::pow(cos(M_PI2*d/this->_epsilon), 2));
    return true;
  }
  return false;
}

template <typename T, typename S>
SmoothIndicatorTriangle2D<T,S>::SmoothIndicatorTriangle2D(Vector<S,2> center, S radius, S mass, S epsilon, S theta)
{
  this->_center = center;
  this->_radius = radius-.5*epsilon;
  this->_theta = theta;
  this->_epsilon = epsilon;
  this->_mass = mass;
  this->_mofi = 0.5 * this->_mass * pow(this->_radius, 2);
  T smallRad = this->_radius * .5;    //sin(30)
  T halfEdge = this->_radius * std::sqrt(3)/2.; // cos(30)

  _PointA[0] = 0.;
  _PointA[1] = this->_radius;
  _PointB[0] = - halfEdge;
  _PointB[1] = - smallRad;
  _PointC[0] = halfEdge;
  _PointC[1] = - smallRad;

  T invEps = 1./this->_epsilon;

  _ab = _PointB - _PointA;
  _ab.normalize(invEps);
  _ab_d = _ab[1]*_PointA[0] - _ab[0]*_PointA[1];
  _bc = _PointC - _PointB;
  _bc.normalize(invEps);
  _bc_d = _bc[1]*_PointB[0] - _bc[0]*_PointB[1];
  _ca = _PointA - _PointC;
  _ca.normalize(invEps);
  _ca_d = _ca[1]*_PointC[0] - _ca[0]*_PointC[1];

  this->_myMin[0] = this->_center[0] - this->_radius - 2*this->_epsilon;
  this->_myMin[1] = this->_center[1] - this->_radius - 2*this->_epsilon;
  this->_myMax[0] = this->_center[0] + this->_radius + 2*this->_epsilon;
  this->_myMax[1] = this->_center[1] + this->_radius + 2*this->_epsilon;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorTriangle2D<T,S>::getMin()
{
  this->_myMin[0] = this->_center[0] - this->_radius - 2*this->_epsilon;
  this->_myMin[1] = this->_center[1] - this->_radius - 2*this->_epsilon;
  return this->_myMin;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorTriangle2D<T,S>::getMax()
{
  this->_myMax[0] = this->_center[0] + this->_radius + 2*this->_epsilon;
  this->_myMax[1] = this->_center[1] + this->_radius + 2*this->_epsilon;
  return this->_myMax;
}

// returns true if x is inside the sphere
template <typename T, typename S>
bool SmoothIndicatorTriangle2D<T,S>::operator()(T output[], const S input[])
{

  T xDist = input[0] - this->_center[0];
  T yDist = input[1] - this->_center[1];

  T ct = std::cos(this->_theta);
  T st = std::sin(this->_theta);
  T x = xDist*ct - yDist*st;
  T y = xDist*st + yDist*ct;

  unsigned short area = 0;

  T dist_a = _bc[1]*x-_bc[0]*y - _bc_d;
  T dist_b = _ca[1]*x-_ca[0]*y - _ca_d;
  T dist_c = _ab[1]*x-_ab[0]*y - _ab_d;

  if (dist_c < 0) {
    area = (area | 100);
  }
  if (dist_a < 0) {
    area = (area | 10);
  }
  if (dist_b < 0) {
    area = (area | 1);
  }

  if (area == 111) {
    output[0] = 1;
    return true;
  }

  if (area == 110 && dist_b < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2));
    return true;
  }

  if (area == 101 && dist_a < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_a), 2));
    return true;
  }

  if (area == 11 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 1 && dist_a < 1 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_a), 2)*std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 10 && dist_b < 1 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2)*std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 100 && dist_b < 1 && dist_a < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2)*std::pow(cos(M_PI2*dist_a), 2));
    return true;
  }

  output[0] = 0;
  return false;
}


template <typename T, typename S>
ParticleIndicatorCuboid2D<T,S>::ParticleIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S mass, S epsilon, S theta)
  : _xLength(xLength),_yLength(yLength)
{
  this->_pos = center;
  this->_epsilon = epsilon;
  this->_theta = theta;
  this->_mass = mass;
  this->_radius = .5*(std::sqrt(std::pow(_xLength, 2)+std::pow(_yLength, 2)));
  this->_mofi = 1./12.*this->_mass*(_xLength*_xLength+_yLength*_yLength);
  this->_myMin[0] = - this->_radius - 2.*this->_epsilon;
  this->_myMin[1] = - this->_radius - 2.*this->_epsilon;
  this->_myMax[0] = this->_radius + 2.*this->_epsilon;
  this->_myMax[1] = this->_radius + 2.*this->_epsilon;
}

template <typename T, typename S>
bool ParticleIndicatorCuboid2D<T,S>::operator()(T output[], const S r[])
{
  T xDist = r[0] - this->_pos[0];
  T yDist = r[1] - this->_pos[1];

  T xL2 = _xLength/2.;
  T yL2 = _yLength/2.;

  // counter-clockwise rotation by _theta=-theta around center
  T ct = std::cos(this->_theta);
  T st = std::sin(this->_theta);

  T x = this->_pos[0] + xDist*ct - yDist*st;
  T y = this->_pos[1] + xDist*st + yDist*ct;

  xDist = fabs(x -this-> _pos[0]);
  yDist = fabs(y -this-> _pos[1]);

  if ( xDist <= xL2 && yDist <= yL2) {
    output[0] = 1.;
    return true;
  }
  if ( xDist > xL2 + this->_epsilon || yDist > yL2 + this->_epsilon ) {
    output[0] = 0.;
    return false;
  }
  if ( xDist < xL2 && (yDist <= yL2 + this->_epsilon  && yDist > yL2) ) {
    output[0] = T( std::pow(cos(M_PI2*(yDist - yL2)/this->_epsilon), 2));
    return true;
  }
  if ( yDist < yL2 && (xDist <= xL2 + this->_epsilon  && xDist > xL2) ) {
    output[0] = T( std::pow(cos(M_PI2*(xDist - xL2)/this->_epsilon), 2));
    return true;
  }
  if ( (xDist <= xL2 + this->_epsilon && xDist > xL2) && (yDist <= yL2 + this->_epsilon && yDist > yL2) ) {
    output[0] = T( (std::pow(cos(M_PI2*(xDist - xL2)/this->_epsilon), 2) *
                    std::pow(cos(M_PI2*(yDist - yL2)/this->_epsilon), 2)) );
    return true;
  }
  output[0] = 0.;
  return false;
}

template <typename T, typename S>
ParticleIndicatorCircle2D<T,S>::ParticleIndicatorCircle2D(Vector<S,2> center, S radius, S mass, S epsilon)
{
  this->_radius = radius;
  this->_epsilon = epsilon;
  this->_pos = center;
  this->_myMin[0] = - this->_radius - 2.*this->_epsilon;
  this->_myMin[1] = - this->_radius - 2.*this->_epsilon;
  this->_myMax[0] = this->_radius + 2.*this->_epsilon;
  this->_myMax[1] = this->_radius + 2.*this->_epsilon;
  this->_mass = mass;
  this->_mofi = 0.5 * this->_mass * pow(this->_radius, 2);
}

// returns true if x is inside the sphere
template <typename T, typename S>
bool ParticleIndicatorCircle2D<T,S>::operator()(T output[], const S input[])
{
  double d;   // distance to the figure
  double distToCenter2 = std::pow(this->_pos[0]-input[0], 2) +
                         std::pow(this->_pos[1]-input[1], 2);
  if ( distToCenter2 >= std::pow(this->_radius + this->_epsilon *0.5, 2)) {
    output[0] = T(0);
    return true;
  } else if ( distToCenter2 <= std::pow(this->_radius - this->_epsilon *0.5, 2)) {
    output[0] = T(1);
    return true;
  } else {
    // d is between 0 and _epsilon
    d = std::sqrt(distToCenter2) - this->_radius + this->_epsilon *0.5;
    output[0] = T(std::pow(cos(M_PI2*d/this->_epsilon), 2));
    return true;
  }
  return false;
}

template <typename T, typename S>
ParticleIndicatorTriangle2D<T,S>::ParticleIndicatorTriangle2D(Vector<S,2> center, S radius, S mass, S epsilon, S theta)
{
  this->_pos = center;
  this->_theta = theta;
  this->_epsilon = epsilon;
  this->_mass = mass;
  this->_radius = radius;
  this->_mofi = 0.5 * this->_mass * pow(this->_radius, 2);
  T smallRad = this->_radius * .5;    //sin(30)
  T halfEdge = this->_radius * std::sqrt(3)/2.; // cos(30)

  _PointA[0] = 0.;
  _PointA[1] = this->_radius;
  _PointB[0] = - halfEdge;
  _PointB[1] = - smallRad;
  _PointC[0] = halfEdge;
  _PointC[1] = - smallRad;

  T invEps = 1./this->_epsilon;

  _ab = _PointB - _PointA;
  _ab.normalize(invEps);
  _ab_d = _ab[1]*_PointA[0] - _ab[0]*_PointA[1];
  _bc = _PointC - _PointB;
  _bc.normalize(invEps);
  _bc_d = _bc[1]*_PointB[0] - _bc[0]*_PointB[1];
  _ca = _PointA - _PointC;
  _ca.normalize(invEps);
  _ca_d = _ca[1]*_PointC[0] - _ca[0]*_PointC[1];

  this->_myMin[0] = - this->_radius - 2.*this->_epsilon;
  this->_myMin[1] = - this->_radius - 2.*this->_epsilon;
  this->_myMax[0] = this->_radius + 2.*this->_epsilon;
  this->_myMax[1] = this->_radius + 2.*this->_epsilon;
}

// returns true if x is inside the triangle
template <typename T, typename S>
bool ParticleIndicatorTriangle2D<T,S>::operator()(T output[], const S input[])
{

  T xDist = input[0] - this->_pos[0];
  T yDist = input[1] - this->_pos[1];

  T ct = std::cos(this->_theta);
  T st = std::sin(this->_theta);
  T x = xDist*ct - yDist*st;
  T y = xDist*st + yDist*ct;

  unsigned short area = 0;

  T dist_a = _bc[1]*x-_bc[0]*y - _bc_d;
  T dist_b = _ca[1]*x-_ca[0]*y - _ca_d;
  T dist_c = _ab[1]*x-_ab[0]*y - _ab_d;

  if (dist_c < 0) {
    area = (area | 100);
  }
  if (dist_a < 0) {
    area = (area | 10);
  }
  if (dist_b < 0) {
    area = (area | 1);
  }

  if (area == 111) {
    output[0] = 1;
    return true;
  }

  if (area == 110 && dist_b < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2));
    return true;
  }

  if (area == 101 && dist_a < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_a), 2));
    return true;
  }

  if (area == 11 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 1 && dist_a < 1 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_a), 2)*std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 10 && dist_b < 1 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2)*std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 100 && dist_b < 1 && dist_a < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2)*std::pow(cos(M_PI2*dist_a), 2));
    return true;
  }

  output[0] = 0;
  return false;
}

template <typename T, typename S>
ParticleIndicatorCustom2D<T,S>::ParticleIndicatorCustom2D(LBconverter<T> const& converter,
			IndicatorF3D<T>& ind,
			Vector<T,2> center,
			T rhoP,
            T epsilon,
            T theta,
			T slice)
            : _converter(converter)
  {
	OstreamManager clout(std::cout,"createIndicatorCustom2D");
	this->_pos = center;
    this->_epsilon = epsilon;
    this->_theta = theta;

    // initialize temporary values
    SmoothBlockIndicator3D<T,olb::descriptors::D3Q19Descriptor> smoothBlock(ind, this->_epsilon);
    int _nX = smoothBlock.getBlockData().getNx();
    int _nY = smoothBlock.getBlockData().getNy();
    int tmpNcells = 0;
    if (slice<ind.getMin()[2] || slice>ind.getMax()[2]) {
    	clout << "ERROR: Forbidden value, slice out of bounds" << std::endl;
    	return;
    }

    // create smoothed blockData
    int tmpZ = int(slice/this->_converter.getLatticeL());
    BlockData2D<T,BaseType> block_tmp(_nX, _nY);
    for (int iX=0; iX < _nX; iX++) {
      for (int iY=0; iY < _nY; iY++) {
        block_tmp.get(iX, iY) = smoothBlock.getBlockData().get(iX, iY, tmpZ);
        if(block_tmp.get(iX, iY) > 0) {
          tmpNcells++;
        }
      }
    }
    this->_blockData = block_tmp;
    T invNcells = 1./tmpNcells;

    // calculate mass and centerpoint for rotation
    this->_mass = rhoP * tmpNcells * _converter.physLength()*_converter.physLength();
    this->_center[0] = 0.0;
    this->_center[1] = 0.0;
    for(int iX= 0; iX < _nX; iX++) {
      for(int iY = 0; iY < _nY; iY++) {
        if(this->_blockData.get(iX,iY) > std::numeric_limits<T>::epsilon()) {
          this->_center[0] += (this->_converter.physLength(iX)) * this->_blockData.get(iX, iY) * invNcells;
          this->_center[1] += (this->_converter.physLength(iY)) * this->_blockData.get(iX, iY) * invNcells;
        }
      }
    }
    this->_latticeCenter[0] = _converter.numCells(this->_center[0]);
    this->_latticeCenter[1] = _converter.numCells(this->_center[1]);

    // calculate moment of inertia
    this->_mofi = 0.;
    for(int iX = 0; iX < _nX; ++iX) {
      for(int iY = 0; iY < _nY; ++iY) {
        if(this->_blockData.get(iX,iY) > std::numeric_limits<T>::epsilon()) {
          T dx = std::abs(this->_converter.physLength(iX) - this->_center[0]);
          T dy = std::abs(this->_converter.physLength(iY) - this->_center[1]);
          this->_mofi += (dx*dx+dy*dy);
        }
      }
    }
    this->_mofi += pow(this->_converter.physLength(1), 4)/ 6.0;
    this->_mofi *= this->_mass*invNcells;

    // calculate min and max from circumradius
    T distance = 0.;
    for(int iX = 0; iX < _nX; iX++) {
	  T x = this->_converter.physLength(iX);
      for(int iY = 0; iY < _nY; iY++) {
	    T y = this->_converter.physLength(iY);
        if(this->_blockData.get(iX,iY) > std::numeric_limits<T>::epsilon()) {
    	  T tmpDist = std::sqrt(std::pow(this->_center[0]-x,2)+std::pow(this->_center[1]-y,2));
    	  if (tmpDist > distance) distance = tmpDist;
        }
      }
    }
    this->_myMin[0] = -this->_epsilon - distance;
    this->_myMin[1] = -this->_epsilon - distance;
    this->_myMax[0] = this->_epsilon + distance;
    this->_myMax[1] = this->_epsilon + distance;
  }

template <typename T, typename S>
bool ParticleIndicatorCustom2D<T,S>::operator() (T output[], const S input[]) {
  // Translation
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T ct = std::cos(this->getTheta());
  T st = std::sin(this->getTheta());

  // counter-clockwise rotation by _theta=-theta around (0/0) and movement from rotation center to local center
  int x = this->_converter.numCells(xDist*ct - yDist*st) + this->_latticeCenter[0];
  int y = this->_converter.numCells(xDist*st + yDist*ct) + this->_latticeCenter[1];

  /// Checking if coordinates are inside the BlockData
  if(x >= 0 && x < _blockData.getNx() && y >= 0 && y < _blockData.getNy()) {
    if(this->_blockData.get(x, y) > std::numeric_limits<T>::epsilon()) {
      output[0] = T(this->_blockData.get(x, y));
      return true;
    }
  }
  output[0] = T(0);
  return false;
}

} // namespace olb

#endif
