/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2016 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Benjamin Förster
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

#ifndef SUPER_LATTICE_INTEGRAL_F_3D_H
#define SUPER_LATTICE_INTEGRAL_F_3D_H

#include<vector>

#include "functors/genericF.h"
#include "functors/blockLatticeIntegralF3D.h"
#include "functors/superBaseF3D.h"
#include "functors/indicator/superIndicatorBaseF3D.h"
#include "functors/indicator/indicatorBaseF3D.h"
#include "functors/indicator/superIndicatorF3D.h"
#include "functors/superLatticeLocalF3D.h"
#include "functors/interpolationF3D.h"
#include "core/superLattice3D.h"
#include "core/vector.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry3D.h"



/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////


/// functor that returns the max in each component of all points of a certain material
template <typename T, typename W = T>
class SuperMin3D final : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperMin3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
  bool operator() (W output[], const int input[]);
};

/// functor that returns the max in each component of all points of a certain material
template <typename T, typename W = T>
class SuperMax3D final : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperMax3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
  bool operator() (W output[], const int input[]);
};


/// sums over all cells of a certain material number
template <typename T, typename W = T>
class SuperSum3D final : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperSum3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
  bool operator() (W output[], const int input[]);
};

/// sums over all cells of a SmoothIndicatorSphere
template <typename T, typename W = T>
class SuperSumIndicator3D final : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  ParticleIndicatorF3D<T,T>& _indicator;
public:
  SuperSumIndicator3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, ParticleIndicatorF3D<T,T>& indicator);
  bool operator() (T output[], const int input[]);
};


/// average over all cells of a certain material number
template <typename T, typename W = T>
class SuperAverage3D final : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperAverage3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
  bool operator() (W output[], const int input[]);
};

template <typename T, typename W = T>
class SuperIntegral3D final : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperIntegral3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
  bool operator() (W output[], const int input[]);
};

/// functor that returns the Lp norm over omega of the the euklid norm of the input functor
/**
 * p == 0: inf norm
 * p >= 1: p norm
 */
template <typename T, typename W = T>
class SuperLpNorm3D : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  SuperIndicatorF3D<T>& _indicatorF;
  int _p;
public:
  SuperLpNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, SuperIndicatorF3D<T>& indicatorF, int p);
  SuperLpNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, std::vector<int> materials, int p);
  SuperLpNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material, int p);
  bool operator() (W output[], const int input[]);
};


/// functor that returns the L1 norm over omega of the the euklid norm of the input functor
template <typename T, typename W = T>
class SuperL1Norm3D final : public SuperLpNorm3D<T,W> {
public:
  SuperL1Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, SuperIndicatorF3D<T>& indicatorF);
  SuperL1Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, std::vector<int> materials);
  SuperL1Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
};


/// functor that returns the L2 norm over omega of the the euklid norm of the input functor
template <typename T, typename W = T>
class SuperL2Norm3D final : public SuperLpNorm3D<T,W> {
public:
  SuperL2Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, SuperIndicatorF3D<T>& indicatorF);
  SuperL2Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, std::vector<int> materials);
  SuperL2Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
};


/// functor that returns the Linf norm over omega of the the euklid norm of the input functor
template <typename T, typename W = T>
class SuperLinfNorm3D final : public SuperLpNorm3D<T,W> {
public:
  SuperLinfNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, SuperIndicatorF3D<T>& indicatorF);
  SuperLinfNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, std::vector<int> materials);
  SuperLinfNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
};

///// functor that returns the Linf norm over omega (with given material) of the the euklid norm of the input functor
///// \todo implement LpNorm with p=0 as infNorm
//template <typename T, typename W = T>
//class SuperLinfNorm3D final : public SuperF3D<T,W> {
//private:
//  SuperF3D<T,W>& _f;
//  SuperGeometry3D<T>& _superGeometry;
//  const int _material;
//public:
//  SuperLinfNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, const int material);
//  bool operator() (W output[], const int input[]);
//};



/// functor counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T>
class SuperGeometryFaces3D final : public GenericF<T,int> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
  const LBconverter<T>& _converter;
  std::vector<BlockGeometryFaces3D<T>* > _blockGeometryFaces;
public:
  SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry, const int material,
                       const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
  std::vector<T> operator() (std::vector<int> input)
  {
    std::vector<T> output(this->getTargetDim(), T());
    operator()(&output[0],&input[0]);
    return output;
  };
};


/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDrag3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
  SuperGeometryFaces3D<T> _faces;
  SuperLatticePhysBoundaryForce3D<T,DESCRIPTOR> _pBoundForce;
  SuperSum3D<T,T> _sumF;
  T _factor;
public:
  SuperLatticePhysDrag3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry3D<T>& superGeometry, const int material,
                         const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise phys force acting on a boundary with a given indicator on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDragIndicator3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  ParticleIndicatorSphere3D<T,T>& _indicator;
public:
  SuperLatticePhysDragIndicator3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                  SuperGeometry3D<T>& superGeometry,
                                  ParticleIndicatorSphere3D<T,T>& indicator,
                                  const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/**
 *  functor to get pointwise phys force acting on a boundary with a given material on local lattice
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysCorrDrag3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysCorrDrag3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             SuperGeometry3D<T>& superGeometry, const int material,
                             const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/**
 *  functor to get the surface integral of a vector field where the vector field
 *  is represented by a SuperLatticeF functor and the surface is a plane (or area of a circle)
 */
template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticeFlux3D final : public SuperLatticeF3D<T, DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _sg;
  /// define the plane by a starting point _origin and either a _normal, or two vectors _u and _v
  Vector<T,3> _u;
  Vector<T,3> _v;
  Vector<T,3> _origin;
  Vector<T,3> _normal;
  /// radius of the plane (optional)
  T _rad;
  /// grid length (default=lattice length)
  T _h;
  /// number of points of the (discretized) plane
  int _vox;
  /// list of materials to check on the plane (optional; default=1)
  std::list<int> _mat;
  /// functor for interpolation
  AnalyticalFfromSuperF3D<T> _analyticalF;

  /// initializes the member variables ie. defines all variables concerning the plane
  void init(SuperLatticeF3D<T, DESCRIPTOR>& f);
  /// checks if point physR (and its direct neighbours) have the material numbers of _mat (default: _mat=1) ie. are inside the domain
  bool checkInside(std::vector<T> physR, int iC);
  /// interpolates the quantity at all points of the area and sums it up
  void calculate(Vector<T,3>& flow, int xDir, int yDir, int xStart=0, int yStart=0);

public:
  /// define plane by
  /// two vectors
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
                     std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
                     T radius = T(), T h = T());
  /// define plane by
  /// two vectors and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
                     std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
                     std::list<int> materials, T radius = T(), T h = T());
  /// define plane by
  /// normal
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
                     std::vector<T>& n, std::vector<T> A, T radius = T(), T h = T() );
  /// define plane by
  /// normal and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
                     std::vector<T>& n, std::vector<T> A, std::list<int> materials,
                     T radius = T(), T h = T());
  /// define plane by
  /// circle
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
                     IndicatorCircle3D<T>& circle, T h = T());
  /// define plane by
  /// circle and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
                     IndicatorCircle3D<T>& circle, std::list<int> materials, T h = T());

  /// returns vector with
  /// output[0]=flux, output[1]=size of the area and output[2..4]=flow vector (ie. vector of summed up quantity)
  /// if quantity has dimension one: output[0] (=flux) is replaced by the force
  bool operator() (T output[], const int input[]);

  //  std::string name() {
  //    return "SuperLatticeFlux3D";
  //  }
  void print(std::string regionName = "", std::string fluxSiScaleName = "", std::string meanSiScaleName = "");
};


template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticePhysPressureFlux3D final : public SuperLatticeF3D<T, DESCRIPTOR> {
private:
  SuperLatticePhysPressure3D<T, DESCRIPTOR> _p;
  SuperLatticeFlux3D<T, DESCRIPTOR> _fluxF;
  mutable OstreamManager clout;

public:
  /// define plane by
  /// two vectors
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysPressureFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
                                 T radius = T(), T h = T());
  /// define plane by
  /// two vectors and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysPressureFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
                                 std::list<int> materials, T radius = T(), T h = T());
  /// define plane by
  /// normal
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysPressureFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A, T radius = T(),
                                 T h = T());
  ///define plane by
  ///normal and material list
  ///a radius (default=0 -> radius will be initialized as diameter of the geometry)
  ///and grid length h (default=latticeL)
  SuperLatticePhysPressureFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A,
                                 std::list<int> materials, T radius = T(), T h = T());
  /// define plane by
  /// circle
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysPressureFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 IndicatorCircle3D<T>& circle, T h = T());
  /// define plane by
  /// circle and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysPressureFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 IndicatorCircle3D<T>& circle, std::list<int> materials,
                                 T h = T());

  /// returns vector with
  /// output[0]=flux, output[1]=size of the area and output[2..4]=flow vector (ie. vector of summed up quantity)
  /// if quantity has dimension one: output[0] (=flux) is replaced by the force
  bool operator() (T output[], const int input[]);
  //  std::string name() {
  //    return "SuperLatticePhysPresssureFlux3D";
  //  }

  void print(std::string regionName = "", std::string fluxSiScaleName = "N", std::string meanSiScaleName = "Pa");
};

template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticePhysVelocityFlux3D final : public SuperLatticeF3D<T, DESCRIPTOR> {
private:
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> _vel;
  SuperLatticeFlux3D<T, DESCRIPTOR> _fluxF;
  mutable OstreamManager clout;

public:
  /// define plane by
  /// two vectors u, v
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  /// vector A is origin
  SuperLatticePhysVelocityFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
                                 T radius = T(), T h = T());
  /// define plane by
  /// two vectors u, v and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  /// vector A is origin
  SuperLatticePhysVelocityFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
                                 std::list<int> materials, T radius = T(), T h = T());
  /// define plane by
  /// normal vector n
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  /// vector A is origin
  SuperLatticePhysVelocityFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A, T radius = T(),
                                 T h = T());
  /// define plane by
  /// normal vector n and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  /// vector A is origin
  SuperLatticePhysVelocityFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A, std::list<int> materials,
                                 T radius = T(), T h = T());
  /// define plane by
  /// circle
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysVelocityFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 IndicatorCircle3D<T>& circle, T h = T());
  /// define plane by
  /// circle and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysVelocityFlux3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T>& converter, SuperGeometry3D<T>& sg,
                                 IndicatorCircle3D<T>& circle,
                                 std::list<int> materials, T h = T());

  /// returns vector with
  /// output[0]=flux, output[1]=size of the area and output[2..4]=flow vector (ie. vector of summed up quantity)
  /// if quantity has dimension one: output[0] (=flux) is replaced by the force
  bool operator() (T output[], const int input[]);
  //  std::string name() {
  //    return "SuperLatticePhysVelocityFlux3D";
  //  }

  void print(std::string regionName = "", std::string fluxSiScaleName = "m^3/s", std::string meanSiScaleName = "m/s");
};


} // end namespace olb

#endif
