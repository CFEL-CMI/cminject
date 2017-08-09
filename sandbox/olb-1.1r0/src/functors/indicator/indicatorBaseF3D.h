/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Albert Mink, Mathias J. Krause, Benjamin Förster
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

#ifndef INDICATOR_BASE_F_3D_H
#define INDICATOR_BASE_F_3D_H


#include <vector>

#include "core/vector.h"
#include "functors/analyticalBaseF.h"
#include "functors/genericF.h"


namespace olb {

/** IndicatorF3D is an application from \f$ \Omega \subset R^3 \to \{0,1\} \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename S>
class IndicatorF3D : public GenericF<bool,S> {
protected:
  IndicatorF3D();
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
public:
  virtual Vector<S,3>& getMin();
  virtual Vector<S,3>& getMax();
  /** \returns false or true and pos. distance if there was one found for a given origin and direction.
   * Mind that the default computation is done by a numerical approximation which searches .. [TODO: CYRIL]
   */
  virtual bool distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC=-1);
  /// returns true and the normal if there was one found for an given origin and direction
  /**
   * (mind that the default computation is done by a numerical approximation which searches .. [TODO])
   */
  virtual bool normal(Vector<S,3>& normal, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC=-1);
  ///Rotate vector around axis by angle theta
  virtual bool rotOnAxis(Vector<S,3>& vec_rot, const Vector<S,3>& vec, const Vector<S,3>& axis, S& theta);
  /// Returns true if `point` is inside a cube with corners `_myMin` and `_myMax`
  bool isInsideBox(Vector<S,3> point);

  /// + Operator (Union)
  IndicatorF3D<S>& operator+(IndicatorF3D<S>& rhs);
  /// - Operator (Without)
  IndicatorF3D<S>& operator-(IndicatorF3D<S>& rhs);
  /// * Operator (Intersection)
  IndicatorF3D<S>& operator*(IndicatorF3D<S>& rhs);
};


template <typename S>
class IndicatorIdentity3D : public IndicatorF3D<S> {
protected:
  IndicatorF3D<S>& _f;
public:
  IndicatorIdentity3D(IndicatorF3D<S>& f);
  bool operator() (bool output[], const S input[]) override;
};


/** SmoothIndicatorF3D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename T, typename S>
class SmoothIndicatorF3D : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D();
  ~SmoothIndicatorF3D() {};
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
  Vector<S,3> _center;
  Vector<S,3> _vel;
  Vector<S,3> _acc;
  Vector<S,6> _acc2;
  Vector<S,3>  _theta;
  Vector<S,3>  _omega;
  Vector<S,3>  _alpha;

  S _mass;
  S _mofi; //Moment of Inertia

  S _epsilon;
  S _radius;
public:
  virtual Vector<S,3>& getMin();
  virtual Vector<S,3>& getMax();
  virtual Vector<S,3>& getCenter();
  virtual Vector<S,3>& getVel();
  virtual Vector<S,3>& getAcc();
  virtual Vector<S,6>& getAcc2();
  virtual Vector<S,3>& getTheta();
  virtual Vector<S,3>& getOmega();
  virtual Vector<S,3>& getAlpha();
  virtual S& getMass();
  virtual S& getMofi();
  virtual S getDiam();
  virtual S getRadius();
  virtual void setCenter(S centerX, S centerY, S centerZ);
  virtual void setTheta(S thetaX, S thetaY, S thetaZ);

  SmoothIndicatorF3D<T,S>& operator+(SmoothIndicatorF3D<T,S>& rhs);
};


template <typename T, typename S>
class SmoothIndicatorIdentity3D : public SmoothIndicatorF3D<T,S> {
protected:
  SmoothIndicatorF3D<T,S>& _f;
public:
  SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};

//TODO remove now unnecessary stuff from smoothIndicator
/** ParticleIndicatorF3D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename T, typename S>
class ParticleIndicatorF3D : public AnalyticalF3D<T,S> {
protected:
  ParticleIndicatorF3D();
  ~ParticleIndicatorF3D() {};
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
  Vector<S,3> _pos;
  Vector<S,3> _vel;
  Vector<S,3> _acc;
  Vector<S,3> _acc2;
  Vector<S,3> _theta;
  Vector<S,3> _omega;
  Vector<S,3> _alpha;
  Vector<S,3> _alpha2;
  Vector<S,3> _mofi;  //moment of inertia
  Vector<S,9> _rotMat;  //saved values of rotation matrix
  S _mass;
  S _epsilon;
  S _radius;

public:
  Vector<S,3>& getMin();
  Vector<S,3>& getMax();
  Vector<S,3>& getPos();
  Vector<S,3>& getVel();
  Vector<S,3>& getAcc();
  Vector<S,3>& getAcc2();
  Vector<S,3>& getTheta();
  Vector<S,3>& getOmega();
  Vector<S,3>& getAlpha();
  Vector<S,3>& getAlpha2();
  Vector<S,3>& getMofi();
  Vector<S,9>& getRotationMat();
  S& getMass();
  S& getRadius();
  S getDiam();

  ParticleIndicatorF3D<T,S>& operator+(ParticleIndicatorF3D<T,S>& rhs);
};

template <typename T, typename S>
class ParticleIndicatorIdentity3D : public ParticleIndicatorF3D<T,S> {
protected:
  SmoothIndicatorF3D<T,S>& _f;
public:
  ParticleIndicatorIdentity3D(ParticleIndicatorF3D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};

}

#endif
