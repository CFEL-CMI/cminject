/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2015 Gilles Zahnd, Mathias J. Krause
 *  Marie-Luise Maier
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

#ifndef FRAME_CHANGE_F_3D_HH
#define FRAME_CHANGE_F_3D_HH

#include<cmath>

#include "frameChangeF3D.h"
#include "frameChangeF2D.h"
#include "core/superLattice3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "utilities/vectorHelpers.h"  // for normalize
#include "geometry/superGeometry3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace olb {

template <typename T>
RotatingLinear3D<T>::RotatingLinear3D(std::vector<T> axisPoint_,
                                      std::vector<T> axisDirection_, T w_, T scale_)
  : AnalyticalF3D<T,T>(3), axisPoint(axisPoint_), axisDirection(axisDirection_),
    w(w_), scale(scale_) { }


template <typename T>
bool RotatingLinear3D<T>::operator()(T output[], const T x[])
{
  output[0] = (axisDirection[1]*(x[2]-axisPoint[2]) - axisDirection[2]*(x[1]-axisPoint[1]))*w*scale;
  output[1] = (axisDirection[2]*(x[0]-axisPoint[0]) - axisDirection[0]*(x[2]-axisPoint[2]))*w*scale;
  output[2] = (axisDirection[0]*(x[1]-axisPoint[1]) - axisDirection[1]*(x[0]-axisPoint[0]))*w*scale;
  return true;
}


template <typename T>
RotatingQuadratic1D<T>::RotatingQuadratic1D(std::vector<T> axisPoint_,
    std::vector<T> axisDirection_, T w_, T scale_, T additive_)
  : AnalyticalF3D<T,T>(1), w(w_), scale(scale_), additive(additive_)
{
  axisPoint.resize(3);
  axisDirection.resize(3);
  for (int i = 0; i < 3; ++i) {
    axisPoint[i] = axisPoint_[i];
    axisDirection[i] = axisDirection_[i];
  }
}


template <typename T>
bool RotatingQuadratic1D<T>::operator()(T output[], const T x[])
{
  T radiusSquared =  ((x[1]-axisPoint[1])*(x[1]-axisPoint[1])
                      +(x[2]-axisPoint[2])*(x[2]-axisPoint[2]))*axisDirection[0]
                     + ((x[0]-axisPoint[0])*(x[0]-axisPoint[0])
                        +(x[2]-axisPoint[2])*(x[2]-axisPoint[2]))*axisDirection[1]
                     + ((x[1]-axisPoint[1])*(x[1]-axisPoint[1])
                        +(x[0]-axisPoint[0])*(x[0]-axisPoint[0]))*axisDirection[2];
  output[0] = scale*w*w/2*radiusSquared+additive;
  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
RotatingForceField3D<T,DESCRIPTOR>::RotatingForceField3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice_, SuperGeometry3D<T>& superGeometry_,
 const LBconverter<T>& converter_, std::vector<T> axisPoint_,
 std::vector<T> axisDirection_, T w_, bool centrifugeForceOn_,
 bool coriolisForceOn_)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice_,3), sg(superGeometry_),
    converter(converter_), axisPoint(axisPoint_), axisDirection(axisDirection_),
    w(w_), centrifugeForceOn(centrifugeForceOn_), coriolisForceOn(coriolisForceOn_),
    velocity(sLattice_,converter_), rho(sLattice_)
{ }


template <typename T, template <typename U> class DESCRIPTOR>
bool RotatingForceField3D<T,DESCRIPTOR>::operator()(T output[], const int x[])
{
  std::vector<T> force(3,0);
  std::vector<T> F_centri(3,0);
  std::vector<T> F_coriolis(3,0);

  if ( this->_sLattice.getLoadBalancer().rank(x[0]) == singleton::mpi().getRank() ) {
    // local coords are given, fetch local cell and compute value(s)
    std::vector<T> physR(3,T());
    this->sg.getCuboidGeometry().getPhysR(&(physR[0]),&(x[0]));

    T scalar =  (physR[0]-axisPoint[0])*axisDirection[0]
                +(physR[1]-axisPoint[1])*axisDirection[1]
                +(physR[2]-axisPoint[2])*axisDirection[2];

    if (centrifugeForceOn) {
      F_centri[0] = w*w*(physR[0]-axisPoint[0]-scalar*axisDirection[0]);
      F_centri[1] = w*w*(physR[1]-axisPoint[1]-scalar*axisDirection[1]);
      F_centri[2] = w*w*(physR[2]-axisPoint[2]-scalar*axisDirection[2]);
    }
    if (coriolisForceOn) {
      T _vel[3];
      (velocity)(_vel,x);
      F_coriolis[0] = -2*w*(axisDirection[1]*_vel[2]-axisDirection[2]*_vel[1]);
      F_coriolis[1] = -2*w*(axisDirection[2]*_vel[0]-axisDirection[0]*_vel[2]);
      F_coriolis[2] = -2*w*(axisDirection[0]*_vel[1]-axisDirection[1]*_vel[0]);
    }
    T _rho[1];
    (rho)(_rho,x);
    force[0] = (_rho[0]*(F_coriolis[0]+F_centri[0]))/converter.physMasslessForce();
    force[1] = (_rho[0]*(F_coriolis[1]+F_centri[1]))/converter.physMasslessForce();
    force[2] = (_rho[0]*(F_coriolis[2]+F_centri[2]))/converter.physMasslessForce();
  }
  return true;
}


template <typename T>
CirclePoiseuille3D<T>::CirclePoiseuille3D(std::vector<T> center, std::vector<T> normal,
    T maxVelocity, T radius, T scale)
  : AnalyticalF3D<T,T>(3), _center(center), _normal(util::normalize(normal)),
    _radius(radius), _maxVelocity(maxVelocity), _scale(scale) { }


template <typename T>
CirclePoiseuille3D<T>::CirclePoiseuille3D(SuperGeometry3D<T>& superGeometry,
    int material, T maxVelocity, T scale)
  : AnalyticalF3D<T,T>(3), _maxVelocity(maxVelocity), _scale(scale)
{
  _center = superGeometry.getStatistics().getCenterPhysR(material);
  std::vector<T> normalTmp;
  normalTmp.push_back(superGeometry.getStatistics().computeDiscreteNormal(material)[0]);
  normalTmp.push_back(superGeometry.getStatistics().computeDiscreteNormal(material)[1]);
  normalTmp.push_back(superGeometry.getStatistics().computeDiscreteNormal(material)[2]);
  _normal = util::normalize(normalTmp);

  // div. by 2 since one of the discrete redius directions is always 0 while the two other are not
  _radius = T();
  for (int iD = 0; iD < 3; iD++) {
    _radius += superGeometry.getStatistics().getPhysRadius(material)[iD]/2.;
  }
}


template <typename T>
bool CirclePoiseuille3D<T>::operator()(T output[], const T x[])
{
  output[0] = _scale*_maxVelocity*_normal[0]*(1.-((x[1]-_center[1])*(x[1]-_center[1])+(x[2]-_center[2])*(x[2]-_center[2]))/_radius/_radius);
  output[1] = _scale*_maxVelocity*_normal[1]*(1.-((x[0]-_center[0])*(x[0]-_center[0])+(x[2]-_center[2])*(x[2]-_center[2]))/_radius/_radius);
  output[2] = _scale*_maxVelocity*_normal[2]*(1.-((x[1]-_center[1])*(x[1]-_center[1])+(x[0]-_center[0])*(x[0]-_center[0]))/_radius/_radius);

  return true;
}


template <typename T>
RectanglePoiseuille3D<T>::RectanglePoiseuille3D(std::vector<T>& x0_, std::vector<T>& x1_,
    std::vector<T>& x2_, std::vector<T>& maxVelocity_)
  : AnalyticalF3D<T,T>(3), clout(std::cout,"RectanglePoiseille3D"), x0(x0_), x1(x1_),
    x2(x2_), maxVelocity(maxVelocity_) { }


template <typename T>
RectanglePoiseuille3D<T>::RectanglePoiseuille3D(SuperGeometry3D<T>& superGeometry_,
    int material_, std::vector<T>& maxVelocity_, T offsetX, T offsetY, T offsetZ)
  : AnalyticalF3D<T,T>(3), clout(std::cout, "RectanglePoiseuille3D"), maxVelocity(maxVelocity_)
{
  std::vector<T> min = superGeometry_.getStatistics().getMinPhysR(material_);
  std::vector<T> max = superGeometry_.getStatistics().getMaxPhysR(material_);
  // set three corners of the plane, move by offset to land on the v=0
  // boundary and adapt certain values depending on the orthogonal axis to get
  // different points
  x0 = min;
  x1 = min;
  x2 = min;
  if ( util::nearZero(max[0]-min[0]) ) { // orthogonal to x-axis
    x0[1] -= offsetY;
    x0[2] -= offsetZ;
    x1[1] -= offsetY;
    x1[2] -= offsetZ;
    x2[1] -= offsetY;
    x2[2] -= offsetZ;

    x1[1] = max[1] + offsetY;
    x2[2] = max[2] + offsetZ;
  } else if ( util::nearZero(max[1]-min[1]) ) { // orthogonal to y-axis
    x0[0] -= offsetX;
    x0[2] -= offsetZ;
    x1[0] -= offsetX;
    x1[2] -= offsetZ;
    x2[0] -= offsetX;
    x2[2] -= offsetZ;

    x1[0] = max[0] + offsetX;
    x2[2] = max[2] + offsetZ;
  } else if ( util::nearZero(max[2]-min[2]) ) { // orthogonal to z-axis
    x0[0] -= offsetX;
    x0[1] -= offsetY;
    x1[0] -= offsetX;
    x1[1] -= offsetY;
    x2[0] -= offsetX;
    x2[1] -= offsetY;

    x1[0] = max[0] + offsetX;
    x2[1] = max[1] + offsetY;
  } else {
    clout << "Error: constructor from material number works only for axis-orthogonal planes" << std::endl;
  }
}


template <typename T>
bool RectanglePoiseuille3D<T>::operator()(T output[], const T x[])
{
  std::vector<T> velocity(3,T());

  // create plane normals and supports, (E1,E2) and (E3,E4) are parallel
  std::vector<std::vector<T> > n(4,std::vector<T>(3,0)); // normal vectors
  std::vector<std::vector<T> > s(4,std::vector<T>(3,0)); // support vectors
  for (int iD = 0; iD < 3; iD++) {
    n[0][iD] =  x1[iD] - x0[iD];
    n[1][iD] = -x1[iD] + x0[iD];
    n[2][iD] =  x2[iD] - x0[iD];
    n[3][iD] = -x2[iD] + x0[iD];
    s[0][iD] = x0[iD];
    s[1][iD] = x1[iD];
    s[2][iD] = x0[iD];
    s[3][iD] = x2[iD];
  }
  for (int iP = 0; iP < 4; iP++) {
    n[iP] = util::normalize(n[iP]);
  }

  // determine plane coefficients in the coordinate equation E_i: Ax+By+Cz+D=0
  // given form: (x-s)*n=0 => A=n1, B=n2, C=n3, D=-(s1n1+s2n2+s3n3)
  std::vector<T> A(4,0);
  std::vector<T> B(4,0);
  std::vector<T> C(4,0);
  std::vector<T> D(4,0);

  for (int iP = 0; iP < 4; iP++) {
    A[iP] = n[iP][0];
    B[iP] = n[iP][1];
    C[iP] = n[iP][2];
    D[iP] = -(s[iP][0]*n[iP][0] + s[iP][1]*n[iP][1] + s[iP][2]*n[iP][2]);
  }

  // determine distance to the walls by formula
  std::vector<T> d(4,0);
  T aabbcc(0);
  for (int iP = 0; iP < 4; iP++) {
    aabbcc = A[iP]*A[iP] + B[iP]*B[iP] + C[iP]*C[iP];
    d[iP] = fabs(A[iP]*x[0]+B[iP]*x[1]+C[iP]*x[2]+D[iP])*pow(aabbcc,-0.5);
  }

  // length and width of the rectangle
  std::vector<T> l(2,0);
  l[0] = d[0] + d[1];
  l[1] = d[2] + d[3];

  T positionFactor = 16.*d[0]/l[0]*d[1]/l[0]*d[2]/l[1]*d[3]/l[1]; // between 0 and 1

  output[0] = maxVelocity[0]*positionFactor;
  output[1] = maxVelocity[1]*positionFactor;
  output[2] = maxVelocity[2]*positionFactor;
  return true;
}


template <typename T>
EllipticPoiseuille3D<T>::EllipticPoiseuille3D(std::vector<T> center, T a, T b, T maxVel)
  : AnalyticalF3D<T,T>(3), clout(std::cout, "EllipticPoiseuille3D"), _center(center),
    _a2(a*a), _b2(b*b), _maxVel(maxVel) { }


template <typename T>
bool EllipticPoiseuille3D<T>::operator()(T output[], const T x[])
{
  T cosOmega2 = std::pow(x[0]-_center[0],2.)/(std::pow(x[0]-_center[0],2.)+std::pow(x[1]-_center[1],2.));
  T rad2 = _a2*_b2 /((_b2-_a2)*cosOmega2 + _a2) ;
  T x2 = (std::pow(x[0]-_center[0],2.)+std::pow(x[1]-_center[1],2.))/rad2;

  //  clout << _center[0] << " " << _center[1] << " " << x[0] << " " << x[1] << " " << std::sqrt(rad2) << " " << std::sqrt(std::pow(x[0]-_center[0],2.)+std::pow(x[1]-_center[1],2.)) << " " << x2 << std::endl;

  output[0] = 0.;
  output[1] = 0.;
  output[2] = _maxVel*(-x2+1);
  if ( util::nearZero(_center[0]-x[0]) && util::nearZero(_center[1]-x[1]) ) {
    output[2] = _maxVel;
  }
  return true;
}


////////////////////////// Coordinate Transformation ////////////////


// constructor defines _referenceVector to obtain angle between 2 vectors,
// vector x has to be turned by angle in mathematical positive sense depending
// to _orientation to lie over _referenceVector
template<typename T, typename S>
AngleBetweenVectors3D<T, S>::AngleBetweenVectors3D(
  std::vector<T> referenceVector, std::vector<T> orientation)
  : AnalyticalF3D<T, S>(1),
    _referenceVector(referenceVector),
    _orientation(orientation)
{
  this->getName() = "const";
}

// operator returns angle between vectors x and _referenceVector
template<typename T, typename S>
bool AngleBetweenVectors3D<T, S>::operator()(T output[], const S x[])
{
  Vector<T, 3> n_x;
  Vector<T, 3> orientation(_orientation);
  T angle = T(0);

  if ( util::nearZero(x[0]) && util::nearZero(x[1]) && util::nearZero(x[2]) ) {
    // if (abs(x[0]) + abs(x[1]) + abs(x[2]) == T()) {
    output[0] = angle;  // angle = 0
    return true;
  } else {
    //Vector<S, 3> x_tmp(x); // check construction
    n_x[0] = x[0];
    n_x[1] = x[1];
    n_x[2] = x[2];
    n_x.normalize();
  }

  Vector<T, 3> n_ref(_referenceVector);
  n_ref.normalize();
  Vector<T, 3> cross = crossProduct3D(n_x, n_ref);
  T n_dot = n_x * n_ref;


  if ( util::nearZero(cross*orientation) ) {
    // std::cout<< "orientation in same plane as x and refvector" << std::endl;
  }
  // angle = Pi, if n_x, n_ref opposite
  if ( util::nearZero(n_x[0]+n_ref[0]) && util::nearZero(n_x[1]+n_ref[1]) && util::nearZero(n_x[2]+n_ref[2]) ) {
    angle = acos(-1);
  }
  // angle = 0, if n_x, n_ref equal
  else if ( util::nearZero(n_x[0]-n_ref[0]) && util::nearZero(n_x[1]-n_ref[1]) && util::nearZero(n_x[2]-n_ref[2]) ) {
    angle = T();
  }
  // angle in (0,Pi) or (Pi,2Pi), if n_x, n_ref not opposite or equal
  else {
    Vector<T, 3> n_cross(cross);
    n_cross.normalize();
    T normal = cross.norm();
    Vector<T, 3> n_orient;

    if ( !util::nearZero(orientation.norm()) ) {
      n_orient = orientation;
      n_orient.normalize();
    } else {
      std::cout << "orientation vector does not fit" << std::endl;
    }
    if ((cross * orientation) > T()) {
      angle = 2*M_PI - atan2(normal, n_dot);
    }
    if ((cross * orientation) < T()) {
      angle = atan2(normal, n_dot);
    }
  }
  output[0] = angle;
  return true;
}


// constructor to obtain rotation round an _rotAxisDirection with angle _alpha
// and _origin
template<typename T, typename S>
RotationRoundAxis3D<T, S>::RotationRoundAxis3D(std::vector<T> origin,
    std::vector<T> rotAxisDirection, T alpha)
  : AnalyticalF3D<T, S>(3),
    _origin(origin),
    _rotAxisDirection(rotAxisDirection),
    _alpha(alpha)
{
  this->getName() = "const";
}

// operator returns coordinates of the resulting point after rotation of x
// with _alpha in math positive sense to _rotAxisDirection through _origin
template<typename T, typename S>
bool RotationRoundAxis3D<T, S>::operator()(T output[], const S x[])
{
  Vector<T, 3> n(_rotAxisDirection);
  if ( !util::nearZero(n.norm()) && n.norm() > 0 ) {
    //std::cout<< "Rotation axis: " << _rotAxisDirection[0] << " " << _rotAxisDirection[1] << " " << _rotAxisDirection[2] << std::endl;
    n.normalize();
    // translation
    Vector<T, 3> x_tmp;
    for (int i = 0; i < 3; ++i) {
      x_tmp[i] = x[i] - _origin[i];
    }
    // rotation
    output[0] = (n[0]*n[0]*(1 - cos(_alpha)) + cos(_alpha)) * x_tmp[0]
                + (n[0]*n[1]*(1 - cos(_alpha)) - n[2]*sin(_alpha))*x_tmp[1]
                + (n[0]*n[2]*(1 - cos(_alpha)) + n[1]*sin(_alpha))*x_tmp[2]
                + _origin[0];

    output[1] = (n[1]*n[0]*(1 - cos(_alpha)) + n[2]*sin(_alpha))*x_tmp[0]
                + (n[1]*n[1]*(1 - cos(_alpha)) + cos(_alpha))*x_tmp[1]
                + (n[1]*n[2]*(1 - cos(_alpha)) - n[0]*sin(_alpha))*x_tmp[2]
                + _origin[1];

    output[2] = (n[2]*n[0]*(1 - cos(_alpha)) - n[1]*sin(_alpha))*x_tmp[0]
                + (n[2]*n[1]*(1 - cos(_alpha)) + n[0]*sin(_alpha))*x_tmp[1]
                + (n[2]*n[2]*(1 - cos(_alpha)) + cos(_alpha))*x_tmp[2]
                + _origin[2];
    return true;
  } else {
    //std::cout<< "Axis is null or NaN" <<std::endl;
    for (int j=0; j<3; j++) {
      output[j] = x[j];
    }
    return true;
  }
}


////////// Cylindric & Cartesian ///////////////


// constructor to obtain Cartesian coordinates of cylindrical coordinates,
// the z-axis direction equals the axis direction of the cylinder
// with _cylinderOrigin
template<typename T, typename S>
CylinderToCartesian3D<T, S>::CylinderToCartesian3D(
  std::vector<T> cylinderOrigin)
  : AnalyticalF3D<T, S>(3),
    _cylinderOrigin(cylinderOrigin)
{
  this->getName() = "CylinderToCartesian3D";
}

// operator returns Cartesian coordinates of cylindrical coordinates,
// input is x[0] = radius, x[1] = phi, x[2] = z
template<typename T, typename S>
bool CylinderToCartesian3D<T, S>::operator()(T output[], const S x[])
{
  PolarToCartesian2D<T, S> pol2cart(_cylinderOrigin);
  pol2cart(output,x);
  output[2] = x[2];
  return true;
}


// constructor to obtain cylindrical coordinates of Cartesian coordinates
template<typename T, typename S>
CartesianToCylinder3D<T, S>::CartesianToCylinder3D(
  std::vector<T> cartesianOrigin, T& angle, std::vector<T> orientation)
  : AnalyticalF3D<T, S>(3),
    _cartesianOrigin(cartesianOrigin),
    _orientation(orientation)
{
  // rotation with angle to (0,0,1)
  std::vector<T> origin(3,T());
  T e3[3]= {T(),T(),T()};
  e3[2] = T(1);
  RotationRoundAxis3D<T, S> rotRAxis(origin, orientation, angle);
  T tmp[3] = {T(),T(),T()};
  rotRAxis(tmp, e3);

  std::vector<T> htmp(tmp,tmp+3);
  _axisDirection = htmp;
}

// constructor to obtain cylindrical coordinates of Cartesian coordinates
template<typename T, typename S>
CartesianToCylinder3D<T, S>::CartesianToCylinder3D(
  std::vector<T> cartesianOrigin, std::vector<T> axisDirection,
  std::vector<T> orientation)
  : AnalyticalF3D<T, S>(3),
    _cartesianOrigin(cartesianOrigin),
    _axisDirection(axisDirection),
    _orientation(orientation)
{
  this->getName() = "const";
}

// constructor to obtain cylindrical coordinates of Cartesian coordinates
template<typename T, typename S>
CartesianToCylinder3D<T, S>::CartesianToCylinder3D(
  T cartesianOriginX, T cartesianOriginY, T cartesianOriginZ,
  T axisDirectionX, T axisDirectionY, T axisDirectionZ,
  T orientationX, T orientationY, T orientationZ)
  : AnalyticalF3D<T, S>(3)
{
  _cartesianOrigin.push_back(cartesianOriginX);
  _cartesianOrigin.push_back(cartesianOriginY);
  _cartesianOrigin.push_back(cartesianOriginZ);

  _axisDirection.push_back(axisDirectionX);
  _axisDirection.push_back(axisDirectionY);
  _axisDirection.push_back(axisDirectionZ);

  _orientation.push_back(orientationX);
  _orientation.push_back(orientationY);
  _orientation.push_back(orientationZ);
}

// operator returns cylindrical coordinates, output[0] = radius(>= 0),
// output[1] = phi in [0, 2Pi), output[2] = z
template<typename T, typename S>
bool CartesianToCylinder3D<T, S>::operator()(T output[], const S x[])
{
  T x_rot[3] = {T(),T(),T()};
  x_rot[0] = x[0];
  x_rot[1] = x[1];
  x_rot[2] = x[2];
//  T x_rot[3] = {x[0], x[1], x[2]};
  Vector<T, 3> e3(T(),T(), T(1));
  Vector<T, 3> axisDirection(_axisDirection);
  Vector<T, 3> orientation(_orientation);

  Vector<T, 3> normal = crossProduct3D(axisDirection, e3);
  Vector<T, 3> normalAxisDir(axisDirection);
  normalAxisDir.normalize();

  // if axis has to be rotated
  if (!( util::nearZero(normalAxisDir[0]) && util::nearZero(normalAxisDir[1]) && util::nearZero(normalAxisDir[2]-1) )) {

    if ( !util::nearZero(norm(orientation)) ) {
      normal = orientation; // if _orientation = 0,0,0 -> segm. fault
    }
    AngleBetweenVectors3D<T,S> angle(util::fromVector3(e3), util::fromVector3(normal));
    T tmp[3] = {T(),T(),T()};
    tmp[0] = axisDirection[0];
    tmp[1] = axisDirection[1];
    tmp[2] = axisDirection[2];
    T alpha[1] = {T()};
    angle(alpha, tmp);

    // rotation with angle alpha to (0,0,1)
    RotationRoundAxis3D<T, S> rotRAxis(_cartesianOrigin, util::fromVector3(normal), -alpha[0]);
    rotRAxis(x_rot, x);
  }
  CartesianToPolar2D<T, S> car2pol(_cartesianOrigin, util::fromVector3(e3), _orientation);
  T x_rot_help[2] = {T(),T()};
  x_rot_help[0] = x_rot[0];
  x_rot_help[1] = x_rot[1];

  T output_help[2] = {T(),T()};
  car2pol(output_help, x_rot_help);

  output[0] = output_help[0];
  output[1] = output_help[1];
  output[2] = x_rot[2] - _cartesianOrigin[2];
  return true;  // output[0] = radius, output[1] = phi, output[2] = z
}

// Read only access to _axisDirection
template<typename T, typename S>
std::vector<T> const& CartesianToCylinder3D<T, S>::getAxisDirection() const
{
  return _axisDirection;
}

// Read and write access to _axisDirection
template<typename T, typename S>
std::vector<T>& CartesianToCylinder3D<T, S>::getAxisDirection()
{
  return _axisDirection;
}

// Read only access to _catesianOrigin
template<typename T, typename S>
std::vector<T> const& CartesianToCylinder3D<T, S>::getCartesianOrigin() const
{
  return _cartesianOrigin;
}

// Read and write access to _catesianOrigin
template<typename T, typename S>
std::vector<T>& CartesianToCylinder3D<T, S>::getCartesianOrigin()
{
  return _cartesianOrigin;
}

// Read and write access to _axisDirection using angle
template<typename T, typename S>
void CartesianToCylinder3D<T, S>::setAngle(T angle)
{
  std::vector<T> Null(3,T());
  T e3[3] = {T(),T(),T()};
  e3[2] = T(1);

  RotationRoundAxis3D<T, S> rotRAxis(Null, _orientation, angle);
  T tmp[3] = {T(),T(),T()};
  rotRAxis(tmp,e3);
  for (int i=0; i<3; ++i) {
    _axisDirection[i] = tmp[i];
  }
}


////////// Spherical & Cartesian ///////////////


// constructor to obtain spherical coordinates of Cartesian coordinates
template<typename T, typename S>
SphericalToCartesian3D<T, S>::SphericalToCartesian3D()
//std::vector<T> sphericalOrigin)
  : AnalyticalF3D<T, S>(3)  //, _sphericalOrigin(sphericalOrigin)
{
  this->getName() = "const";
}

// operator returns Cartesian coordinates of spherical coordinates
// (x[0] = radius, x[1] = phi, x[2] = theta)
// phi in x-y-plane, theta in x-z-plane, z axis as orientation
template<typename T, typename S>
bool SphericalToCartesian3D<T, S>::operator()(T output[], const S x[])
{
  output[0] = x[0]*sin(x[2])*cos(x[1]);
  output[1] = x[0]*sin(x[2])*sin(x[1]);
  output[2] = x[0]*cos(x[2]);
  return true;
}


template<typename T, typename S>
CartesianToSpherical3D<T, S>::CartesianToSpherical3D(
  std::vector<T> cartesianOrigin, std::vector<T> axisDirection)
  : AnalyticalF3D<T, S>(3),
    _cartesianOrigin(cartesianOrigin),
    _axisDirection(axisDirection)
{
  this->getName() = "const";
}

// operator returns spherical coordinates output[0] = radius(>= 0),
// output[1] = phi in [0, 2Pi), output[2] = theta in [0, Pi]
template<typename T, typename S>
bool CartesianToSpherical3D<T, S>::operator()(T output[], const S x[])
{
  T x_rot[3] = {T(),T(),T()};
  x_rot[0] = x[0];
  x_rot[1] = x[1];
  x_rot[2] = x[2];

  Vector<T,3> axisDirection(_axisDirection);
  Vector<T,3> e3(T(), T(), T(1));
  Vector<T,3> normal = crossProduct3D(axisDirection,e3);

  Vector<T,3> normalAxisDir = axisDirection;
  normalAxisDir.normalize();
  Vector<T,3> cross;
  // if axis has to be rotated
  if ( !( util::nearZero(normalAxisDir[0]) && util::nearZero(normalAxisDir[1]) && util::nearZero(normalAxisDir[2]-1) ) ) {
    AngleBetweenVectors3D<T, S> angle(util::fromVector3(e3), util::fromVector3(normal));
    T tmp[3] = {T(),T(),T()};
    for (int i=0; i<3; ++i) {
      tmp[i] = axisDirection[i];
    }
    T alpha[1] = {T(0)};
    angle(alpha,tmp);
    // cross is rotation axis
    cross = crossProduct3D(e3, axisDirection);
    // rotation with angle alpha to (0,0,1)
    RotationRoundAxis3D<T, S> rotRAxis(_cartesianOrigin, util::fromVector3(normal), -alpha[0]);
    rotRAxis(x_rot,x);
  }

  CartesianToPolar2D<T, S> car2pol(_cartesianOrigin, util::fromVector3(e3), util::fromVector3(cross));
  //  y[0] = car2pol(x_rot)[0];
  Vector<T, 3> difference;

  for (int i=0; i<3; ++i) {
    difference[i] = x_rot[i] - _cartesianOrigin[i];
  }

  T distance = T();
  for (int i = 0; i <= 2; ++i) {
    distance += difference[i]*difference[i];
  }
  distance = sqrt(distance);

  car2pol(output, x_rot);
  output[0] = distance;
  output[2] = acos(x[2]/output[0]);  // angle between z-axis and r-vector

  return true;  // output[0] = radius, output[1] = phi, output[2] = theta
}


////////// Magnetic Field ///////////////


// Magnetic field that creates magnetization in wire has to be orthogonal to
// the wire. To calculate the magnetic force on particles around a cylinder
// (J. Lindner et al. / Computers and Chemical Engineering 54 (2013) 111-121)
// Fm = mu0*4/3.*PI*radParticle^3*_Mp*radWire^2/r^3 *
//   [radWire^2/r^2+cos(2*theta), sin(2*theta), 0]
template<typename T, typename S>
MagneticForceFromCylinder3D<T, S>::MagneticForceFromCylinder3D(
  CartesianToCylinder3D<T, S>& car2cyl, T length, T radWire, T cutoff, T Mw, T Mp,
  T radParticle, T mu0)
  : AnalyticalF3D<T, S>(3),
    _car2cyl(car2cyl),
    _length(length),
    _radWire(radWire),
    _cutoff(cutoff),
    _Mw(Mw),
    _Mp(Mp),
    _radParticle(radParticle),
    _mu0(mu0)
{
  // Field direction H_0 parallel to fluid flow direction, if not: *(-1)
  _factor = -_mu0*4/3*M_PI*pow(_radParticle, 3)*_Mp*_Mw*pow(_radWire, 2);
}

template<typename T, typename S>
bool MagneticForceFromCylinder3D<T, S>::operator()(T output[], const S x[])
{

  Vector<T, 3> magneticForcePolar;
  T outputTmp[3] = {T(), T(), T()};
  Vector<T, 3> normalAxis(_car2cyl.getAxisDirection());
  normalAxis.normalize();

  Vector<T, 3> relPosition;
  relPosition[0] = (x[0] - _car2cyl.getCartesianOrigin()[0]);
  relPosition[1] = (x[1] - _car2cyl.getCartesianOrigin()[1]);
  relPosition[2] = (x[2] - _car2cyl.getCartesianOrigin()[2]);

  T tmp[3] = {T(),T(),T()};
  _car2cyl(tmp,x);
  T rad = tmp[0];
  T phi = tmp[1];

  T test = relPosition * normalAxis;

  if ( (rad > _radWire && rad <= _cutoff*_radWire) &&
       (T(0) <= test && test <= _length) ) {

    magneticForcePolar[0] = _factor/pow(rad, 3)
                            *(pow(_radWire, 2)/pow(rad, 2) + cos(2*phi));
    magneticForcePolar[1] = _factor/pow(rad, 3)*sin(2*phi);

    // changes radial and angular force components to cartesian components.
    outputTmp[0] = magneticForcePolar[0]*cos(phi)
                   - magneticForcePolar[1]*sin(phi);
    outputTmp[1] = magneticForcePolar[0]*sin(phi)
                   + magneticForcePolar[1]*cos(phi);

    // if not in standard axis direction
    if ( !( util::nearZero(normalAxis[0]) && util::nearZero(normalAxis[1]) && util::nearZero(normalAxis[2]-1) ) ) {
      Vector<T, 3> e3(T(), T(), T(1));
      Vector<T, 3> orientation =
        crossProduct3D(Vector<T,3>(_car2cyl.getAxisDirection()),e3);

      AngleBetweenVectors3D<T,S> angle(util::fromVector3(e3), util::fromVector3(orientation));
      T alpha[1] = {T()};
      T tmp2[3] = {T(),T(),T()};
      for (int i = 0; i<3; ++i) {
        tmp2[i] = _car2cyl.getAxisDirection()[i];
      }
      angle(alpha,tmp2);

      std::vector<T> origin(3, T());
      RotationRoundAxis3D<T, S> rotRAxis(origin, util::fromVector3(orientation), alpha[0]);
      rotRAxis(output,outputTmp);
    } else {
      output[0] = outputTmp[0];
      output[1] = outputTmp[1];
      output[2] = T();
    }
  } else {
    output[0] = outputTmp[0];
    output[1] = outputTmp[1];
    output[2] = outputTmp[1];
  }
  return true;
}


} // end namespace olb
#endif
