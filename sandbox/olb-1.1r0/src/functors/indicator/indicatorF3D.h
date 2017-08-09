/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Albert Mink
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

#ifndef INDICATOR_F_3D_H
#define INDICATOR_F_3D_H

#include<vector>
#include "indicatorBaseF3D.h"
#include "io/xmlReader.h"

#include "core/blockData3D.h"
#include "core/units.h"

/** \file
 * This file contains indicator functions. These return 1 if the given
 * coordinates are inside, and 0 if they are outside of the defined set.
 * Implemented are :
 - Sphere
 - Cylinder
 - Cone
 - Pipe (not yet)
 - Cube (not yet)
 - Cuboid
 - Circle

 * The smoothIndicator functors return values in [0,1]. In particular there is
 * an epsilon enclosure of the set, wherein the return values are smooth and do
 * not jump from 0 to 1.

 Boolean operators allow to create unions and intersections. They can be used
 for example for initialization of a SuperGeometry.
*/

namespace olb {

template<typename T> class IndicatorF3D;
template<typename T, typename S> class SmoothIndicatorF3D;


/// indicator function for a 3D circle
template <typename S>
class IndicatorCircle3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center;
  Vector<S,3> _normal;
  S _radius2;
public:
  IndicatorCircle3D(Vector<S,3> center, Vector<S,3> normal, S radius);
  IndicatorCircle3D(S center0, S center1, S center2, S normal0, S normal1,
                    S normal2, S radius);
  bool operator() (bool output[], const S input[]);
  Vector<S,3> const& getCenter() const;
  Vector<S,3> const& getNormal() const;
  S getRadius() const;
  //virtual bool distance(S& distance, Vector<S,3> origin, Vector<S,3> direction, int iC=-1);
};



/// indicator function for a 3D-sphere
template <typename S>
class IndicatorSphere3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center;
  S _radius2;
public:
  IndicatorSphere3D(Vector<S,3> center, S radius);
  IndicatorSphere3D(const IndicatorSphere3D&);
  bool operator() (bool output[], const S input[]) override;
  virtual bool distance(S& distance, const Vector<S,3>& origin,
                        const Vector<S,3>& direction, int iC=-1) override;
};

/// indicator function for a layer
template <typename S>
class IndicatorLayer3D : public IndicatorF3D<S> {
private:
  IndicatorF3D<S>& _indicatorF;
  S _layerSize;
public:
  IndicatorLayer3D(IndicatorF3D<S>& indicatorF, S layerSize);
  bool operator() (bool output[], const S input[]);
};

/// indicator function for a 3d-cylinder
template <typename S>
class IndicatorCylinder3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius2;
  void init();
public:
  IndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> center2, S radius);
  // TODO: eps??
  IndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> normal, S radius, S eps);
  IndicatorCylinder3D(IndicatorCircle3D<S> const& circleF, S eps);
  bool operator() (bool output[], const S input[]);
};

/// indicator function for a 3d frustum
template <typename S>
class IndicatorCone3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius1;
  S _radius2; // The 2nd radius is optional: if not defined, _center2 is the vertex of the cone
public:
  IndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2, S radius1, S radius2=0);
  bool operator() (bool output[], const S input[]);
};

/** indicator function for a 3d-cuboid, parallel to the planes x=0, y=0, z=0.
 * \param extend must have only positive elements
 * \param xLength must be positive
 */
template <typename S>
class IndicatorCuboid3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center;
  S _xLength;
  S _yLength;
  S _zLength;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboid3D(Vector<S,3> extend, Vector<S,3> origin);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboid3D(S xlength, S ylength, S zlength, Vector<S,3> center);
  /// returns true if input is inside, otherwise false
  bool operator() (bool output[], const S input[]);
};

template <typename S>
class IndicatorCuboidOLD3D : public IndicatorF3D<S> {
private:
  std::vector<S> _center;
  S _xLength;
  S _yLength;
  S _zLength;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboidOLD3D(std::vector<S> extend, std::vector<S> origin);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboidOLD3D(S xlength, S ylength, S zlength, std::vector<S> center);
  /// returns true if input is inside, otherwise false
  bool operator() (bool output[], const S input[]);
};


/// indicator function for a 3d-parallelepiped (including any cuboid)
template <typename S>
class IndicatorParallelepiped3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _origin;   // A vertex of the parallelepiped
  Vector<S,3> _corner1;
  Vector<S,3> _corner2;
  Vector<S,3> _corner3;  // Those three vertex have a common edge with the origin.
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _normI;
  S _normJ;
  S _normK;
public:
  IndicatorParallelepiped3D(Vector<S,3> origin, Vector<S,3> corner1,
                            Vector<S,3> corner2, Vector<S,3> corner3);
  bool operator() (bool output[], const S input[]);
};

/////////creatorFunctions//////////////////////


// creator function for geometric primitives
template <typename S>
IndicatorCircle3D<S>* createIndicatorCircle3D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorSphere3D<S>* createIndicatorSphere3D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorCylinder3D<S>* createIndicatorCylinder3D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorCone3D<S>* createIndicatorCone3D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorCuboid3D<S>* createIndicatorCuboid3D(XMLreader const& params, bool verbose=false);

// arithmetic creator functions
template <typename S>
IndicatorF3D<S>* createIndicatorUnion3D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorF3D<S>* createIndicatorWithout3D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorF3D<S>* createIndicatorIntersection3D(XMLreader const&params, bool verbose=false);

// godfather
template <typename S>
IndicatorF3D<S>* createIndicatorF3D(XMLreader const& params, bool verbose=false);



///////////////////////////SmoothIndicatorF/////////////////////////////////////

/// implements a smooth sphere in 3D with an _epsilon sector
template<typename T, typename S>
class SmoothIndicatorSphere3D: public SmoothIndicatorF3D<T, S> {
private:
  S _innerRad;
  S _outerRad;
  S _epsilon;
public:
  SmoothIndicatorSphere3D(Vector<S, 3> center, S radius, S epsilon, S mass);
  SmoothIndicatorSphere3D(const SmoothIndicatorSphere3D<T, S>& rhs);
  bool operator()(T output[], const S input[]);
  Vector<S, 3>& getCenter();
  S getRadius();
  S getDiam();
};

/// implements a smooth cylinder in 3D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCylinder3D : public SmoothIndicatorF3D<T,S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius2;
  S _epsilon;
public:
  SmoothIndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> center2,
                            S radius, S epsilon);
  bool operator() (T output[], const S input[]);
};

/// implements a smooth cone in 3D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCone3D : public SmoothIndicatorF3D<T,S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius1;
  S _radius2; /**< The 2nd radius is optional: if not defined, _center2 is the vertex of the cone */
  S _epsilon;
public:
  SmoothIndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2,
                        S radius1, S radius2, S epsilon);
  bool operator() (T output[], const S input[]);
};



///////////////////////////ParticleIndicatorF/////////////////////////////////////

/// implements a smooth sphere in 3D with an _epsilon sector for particle simulations
template<typename T, typename S>
class ParticleIndicatorSphere3D: public ParticleIndicatorF3D<T, S> {
public:
  ParticleIndicatorSphere3D(Vector<S, 3> center, S radius, S epsilon, S mass);
  bool operator()(T output[], const S input[]);
};

/** implements a smooth particle cuboid in 3D with an _epsilon sector.
 * TODO construct by density
 * TODO rotation seems weird
 */
template <typename T, typename S>
class ParticleIndicatorCuboid3D : public ParticleIndicatorF3D<T, S> {
private:
  S _xLength;
  S _yLength;
  S _zLength;
public:
  ParticleIndicatorCuboid3D(Vector<S,3> center, S xLength, S yLength, S zLength, S mass, S epsilon, Vector<S,3> theta);
  bool operator()(T output[],const S x[]);
};

/** implements a smooth particle of shape given by in indicator (e.g. STL) in 3D with an _epsilon sector.
 * TODO construct by density
 * TODO check correctness of center and mofi
 */
template<typename T, typename S>
class ParticleIndicatorCustom3D : public ParticleIndicatorF3D<T, S> {
private:
  // _center is the local center, _startPos the center at the start
  Vector<T,3> _center;
  // _latticeCenter gives the center in local lattice coordinates
  Vector<int,3> _latticeCenter;
  BlockData3D<T, T> _blockData;
  LBconverter<T> const& _converter;

public:
  // for now epsilon has to be chosen to be latticeL
  ParticleIndicatorCustom3D(LBconverter<T> const& converter, IndicatorF3D<T>& ind, Vector<T,3> center, T rhoP, T epsilon, Vector<T,3> theta);
  bool operator() (T output[], const S input[]);
};

}

#endif

