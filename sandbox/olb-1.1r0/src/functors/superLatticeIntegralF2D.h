/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause,
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

#ifndef SUPER_LATTICE_INTEGRAL_F_2D_H
#define SUPER_LATTICE_INTEGRAL_F_2D_H

#include<vector>
#include<cmath>

#include "functors/genericF.h"
#include "functors/superBaseF2D.h"
#include "functors/superLatticeLocalF2D.h"
#include "functors/interpolationF2D.h"
#include "core/superLattice2D.h"
#include "core/vector.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////


/// functor that returns the max in each component of all points of a certain material
template <typename T>
class SuperMax2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperMax2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
             const int material);
  bool operator() (T output[], const int input[]) override;
};


/// functor that returns the min in each component of all points of a certain material
template <typename T>
class SuperMin2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperMin2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
             const int material);
  bool operator() (T output[], const int input[]) override;
};

/// sums over all cells of a certain material number
template <typename T>
class SuperSum2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperSum2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
             const int material);
  bool operator() (T output[], const int input[]) override;
};


/// sums over all cells of a certain indicator
template <typename T>
class SuperSumIndicator2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  ParticleIndicatorF2D<T,T>& _indicator;
public:
  SuperSumIndicator2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
                      ParticleIndicatorF2D<T,T>& indicator);
  bool operator() (T output[], const int input[]) override;
};


template <typename T>
class SuperIntegral2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperIntegral2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
                  const int material);
  bool operator() (T output[], const int input[]) override;
};


/// functor that returns componentwise the l1 norm
template <typename T>
class SuperL1Norm2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperL1Norm2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
                const int material);
  bool operator() (T output[], const int input[]) override;
};

/// functor that returns the L2 norm over omega (with given material) of the the euklid norm of the input functor
template <typename T>
class SuperL2Norm2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperL2Norm2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
                const int material);
  bool operator() (T output[], const int input[]) override;
};


/// functor that returns the Linf norm over omega (with given material) of the the euklid norm of the input functor
template <typename T>
class SuperLinfNorm2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLinfNorm2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
                  const int material);
  bool operator() (T output[], const int input[]) override;
};

/// functor that returns componentwise the squared l2-norm
template <typename T>
class SuperL222D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperL222D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
             const int material);
  bool operator() (T output[], const int input[]) override;
};


/// functor counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T>
class SuperGeometryFaces2D final : public GenericF<T,int> {
private:
  SuperGeometry2D<T>&   _superGeometry;
  const int             _material;
  const LBconverter<T>& _converter;
public:
  SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry, const int material,
                       const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/// functor counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T>
class SuperGeometryFacesIndicator2D final : public GenericF<T,int> {
private:
  SuperGeometry2D<T>&   _superGeometry;
  SmoothIndicatorCircle2D<T,T>& _indicator;
  const int             _material;
  const LBconverter<T>& _converter;
public:
  SuperGeometryFacesIndicator2D(SuperGeometry2D<T>& superGeometry, SmoothIndicatorCircle2D<T,T>& indicator, const int material,
                                const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDrag2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysDrag2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry2D<T>& superGeometry, const int material,
                         const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise phys force acting on a boundary with a given indicator on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDragIndicator2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  ParticleIndicatorF2D<T,T>& _indicator;
public:
  SuperLatticePhysDragIndicator2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                  SuperGeometry2D<T>& superGeometry,
                                  ParticleIndicatorF2D<T,T>& indicator,
                                  const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise phys force acting on a boundary with a given indicator on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDragIndicator2D_2 final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  SmoothIndicatorF2D<T,T>& _indicator;
public:
  SuperLatticePhysDragIndicator2D_2(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                    SuperGeometry2D<T>& superGeometry,
                                    SmoothIndicatorF2D<T,T>& indicator,
                                    const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/**
 *  functor to get pointwise phys force acting on a boundary with a given material on local lattice
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysCorrDrag2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysCorrDrag2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                             SuperGeometry2D<T>& superGeometry, const int material,
                             const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};

/**
 *  functor to get the surface integral of a vector field, where the vector field
 *  is represented by a SuperLatticeF functor and the surface is a line
 */
template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticeFlux2D final : public SuperLatticeF2D<T, DESCRIPTOR> {
protected:
  SuperGeometry2D<T>& _sg;
  /// define the line by a vectors u or normal and point A
  Vector<T,2> _u, _origin, _normal;
  /// rad is radius of the line, h is the grid length
  T _rad, _h;
  /// number of points of the (discretized) line
  int _vox;
  /// list of materials
  std::list<int> _mat;
  /// functor for interpolation
  AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> _analyticalF;

  /// initializes the member variables ie. defines all variables concerning the line
  void init(SuperLatticeF2D<T, DESCRIPTOR>& f);
  /// checks if point physR (and its direct neighbours) have the material numbers of _mat (default: _mat=1) ie. are inside the domain
  bool checkInside(std::vector<T> p, int iC);
  /// interpolates the quantity at all points of the area and sums it up
  void calculate(std::vector<T>& flow, int xStart, int xDir=0);

public:
  /// define line by
  /// a vectors
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux2D(SuperLatticeF2D<T, DESCRIPTOR>& f, SuperGeometry2D<T>& sg,
                     std::vector<T>& n, std::vector<T> A, T radius = T(), T h = T());
  /// define line by
  /// normal
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticeFlux2D(SuperLatticeF2D<T, DESCRIPTOR>& f, SuperGeometry2D<T>& sg,
                     std::vector<T>& n, std::vector<T> A, std::list<int> materials,
                     T radius = T(), T h = T());

  /// returns vector with
  /// output[0]=flux, output[1]=size of the area and output[2..3]=flow vector (ie. vector of summed up quantity)
  /// if quantity has dimension one: output[0] (=flux) is replaced by the force
  bool operator() (T output[], const int input[]);

  std::string name()
  {
    return "SuperLatticeFlux2D";
  }

  void print(std::string regionName = "", std::string fluxSiScaleName = "", std::string meanSiScaleName = "");
};




template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticePhysPressureFlux2D final : public SuperLatticeF2D<T, DESCRIPTOR> {
private:
  SuperLatticePhysPressure2D<T, DESCRIPTOR> _p;
  SuperLatticeFlux2D<T, DESCRIPTOR> _fluxF;
  mutable OstreamManager clout;

public:
  /// define line by
  /// a normal
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysPressureFlux2D(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T> const& converter, SuperGeometry2D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A, T radius = T(),
                                 T h = T());
  /// define line by
  /// normal and material list
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// grid length h (default=latticeL)
  SuperLatticePhysPressureFlux2D(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T> const& converter, SuperGeometry2D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A, std::list<int> materials,
                                 T radius = T(), T h = T());


  /// returns vector with
  /// output[0]=flux, output[1]=size of the area and output[2..3]=flow vector (ie. vector of summed up quantity)
  /// if quantity has dimension one: output[0] (=flux) is replaced by the force
  bool operator() (T output[], const int input[]);

  std::string name()
  {
    return "SuperLatticePhysPressureFlux2D";
  }

  void print(std::string regionName = "", std::string fluxSiScaleName = "N", std::string meanSiScaleName = "Pa");
};




template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticePhysVelocityFlux2D final : public SuperLatticeF2D<T, DESCRIPTOR> {
private:
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> _vel;
  SuperLatticeFlux2D<T, DESCRIPTOR> _fluxF;
  mutable OstreamManager clout;

public:
  /// define line by
  /// a vectors
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysVelocityFlux2D(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T> const& converter, SuperGeometry2D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A, T radius = T(),
                                 T h = T());
  /// normal
  /// a radius (default=0 -> radius will be initialized as diameter of the geometry)
  /// and grid length h (default=latticeL)
  SuperLatticePhysVelocityFlux2D(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                                 LBconverter<T> const& converter, SuperGeometry2D<T>& sg,
                                 std::vector<T>& n, std::vector<T> A, std::list<int> materials,
                                 T radius = T(), T h = T());

  /// returns vector with
  /// output[0]=flux, output[1]=size of the area and output[2..3]=flow vector (ie. vector of summed up quantity)
  /// if quantity has dimension one: output[0] (=flux) is replaced by the force
  bool operator() (T output[], const int input[]);

  std::string name()
  {
    return "SuperLatticePhysVelocityFlux2D";
  }

  void print(std::string regionName = "", std::string fluxSiScaleName = "m^2/s", std::string meanSiScaleName = "m/s");
};

} // end namespace olb

#endif
