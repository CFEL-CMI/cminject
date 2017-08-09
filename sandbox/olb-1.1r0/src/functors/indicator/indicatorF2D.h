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

#ifndef INDICATOR_F_2D_H
#define INDICATOR_F_2D_H

#include<vector>
#include "indicatorBaseF2D.h"
#include "io/xmlReader.h"

#include "core/blockData2D.h"
#include "core/units.h"
#include "indicatorBaseF3D.h"

/** \file
 * This file contains indicator functions. These return 1 if the given
 * coordinates are inside, and 0 if they are outside of the defined set.
 * Implemented are :
 - Cuboid
 - Circle

 * The smoothIndicator functors return values in [0,1]. In particular there is
 * an epsilon enclosure of the set, wherein the return values are smooth and do
 * not jump from 0 to 1.

 Boolean operators allow to create unions and intersections. They can be used
 for example for initialization of a SuperGeometry.
*/

namespace olb {



/// indicator function for a 2D-cuboid, parallel to the planes x=0, y=0;
/// theta rotates cuboid around its center, theta in radian measure
template <typename S>
class IndicatorCuboid2D : public IndicatorF2D<S> {
private:
  Vector<S,2> _center;
  S _xLength;
  S _yLength;
  S _theta;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboid2D(Vector<S,2> extend, Vector<S,2> origin, S theta=0);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboid2D(S xlength, S ylength, Vector<S,2> center=S(), S theta=0);
  /// returns true if input is inside, otherwise false
  bool operator() (bool output[], const S input[]);
};


/// indicator function for a 2D circle
template <typename S>
class IndicatorCircle2D : public IndicatorF2D<S> {
private:
  Vector<S,2> _center;
  S _radius2;
public:
  IndicatorCircle2D(Vector<S,2> center, S radius);
  bool operator() (bool output[], const S input[]) override;
  virtual bool distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction,  int iC=-1) override;
  virtual bool normal(Vector<S,2>& normal, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC=-1) override;
};


/////////creatorFunctions//////////////////////
template <typename S>
IndicatorCuboid2D<S>* createIndicatorCuboid2D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorCircle2D<S>* createIndicatorCircle2D(XMLreader const& params, bool verbose=false);


///////////////////////////SmoothIndicatorF/////////////////////////////////////

/** implements a smooth cuboid in 2D with an _epsilon sector.
 * \param mass    TODO
 * \param epsilon
 * \param theta   TODO
 */
template <typename T, typename S>
class SmoothIndicatorCuboid2D : public SmoothIndicatorF2D<T,S> {
private:
  S _xLength;
  S _yLength;
public:
  SmoothIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S mass, S epsilon, S theta=0);
  bool operator()(T output[],const S x[]);
  /// \return radius of a circle at center, that contains the object
  S getRadius();
  S getDiam();
  Vector<S,2>& getMin();
  Vector<S,2>& getMax();
};


/// implements a smooth circle in 2D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCircle2D : public SmoothIndicatorF2D<T,S> {
public:
  SmoothIndicatorCircle2D(Vector<S,2> center, S radius, S mass, S epsilon);
  bool operator() (T output[], const S input[]);
  Vector<S,2>& getMin();
  Vector<S,2>& getMax();
};

/// implements a smooth triangle in 2D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorTriangle2D : public SmoothIndicatorF2D<T, S> {
private:
  /// Eckpunkte des Dreiecks
  Vector<S, 2> _PointA, _PointB, _PointC;
  /// Verbindungsvektoren _ab von _A nach _B, etc.
  Vector<S, 2> _ab, _bc, _ca;
  /// normal on _ab * _A  = _ab_d
  S _ab_d, _bc_d, _ca_d;

public:
  SmoothIndicatorTriangle2D(Vector<S,2> center, S radius, S mass, S epsilon, S theta);
  bool operator() (T output[], const S input[]);
  Vector<S,2>& getMin();
  Vector<S,2>& getMax();

};


///////////////////////////ParticleIndicatorF/////////////////////////////////////

/** implements a smooth particle cuboid in 2D with an _epsilon sector.
 * TODO construct by density
 * TODO rotation seems weird
 */
template <typename T, typename S>
class ParticleIndicatorCuboid2D : public ParticleIndicatorF2D<T,S> {
private:
  S _xLength;
  S _yLength;
public:
  ParticleIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S mass, S epsilon, S theta=0);
  bool operator()(T output[],const S x[]);
};

/** implements a smooth particle circle in 2D with an _epsilon sector.
 */
template <typename T, typename S>
class ParticleIndicatorCircle2D : public ParticleIndicatorF2D<T,S> {
public:
  ParticleIndicatorCircle2D(Vector<S,2> center, S radius, S mass, S epsilon);
  bool operator() (T output[], const S input[]);
};

/** implements a smooth particle triangle in 2D with an _epsilon sector, constructed from circumradius
 * TODO correct mofi computation -> 2D mofi required
 * TODO generic constructor with angles
 */
template <typename T, typename S>
class ParticleIndicatorTriangle2D : public ParticleIndicatorF2D<T, S> {
private:
  /// Eckpunkte des Dreiecks
  Vector<S, 2> _PointA, _PointB, _PointC;
  /// Verbindungsvektoren _ab von _A nach _B, etc.
  Vector<S, 2> _ab, _bc, _ca;
  /// normal on _ab * _A  = _ab_d
  S _ab_d, _bc_d, _ca_d;

public:
  ParticleIndicatorTriangle2D(Vector<S,2> center, S radius, S mass, S epsilon, S theta);
  bool operator() (T output[], const S input[]);
};

template<typename T, typename S>
class ParticleIndicatorCustom2D : public ParticleIndicatorF2D<T, S> {
private:
  // _center is the local center, _startPos the center at the start
  Vector<T,2> _center;
  // _latticeCenter gives the center in local lattice coordinates
  Vector<int,2> _latticeCenter;
  BlockData2D<T, T> _blockData;
  LBconverter<T> const& _converter;

public:
  ParticleIndicatorCustom2D(LBconverter<T> const& converter, IndicatorF3D<T>& ind, Vector<T,2> center, T rhoP, T epsilon, T theta, T slice);
  bool operator() (T output[], const S input[]);
};

}

#endif

