/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Albert Mink, Mathias J. Krause
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

/** \file
 * The description of a blockLatticeReduction -- header file.
 */

#ifndef BLOCK_LATTICE_REDUCTION_2D_H
#define BLOCK_LATTICE_REDUCTION_2D_H

#include "io/ostreamManager.h"
#include "core/blockData2D.h"
#include "core/vector.h"
#include "functors/blockBaseF2D.h"
#include "functors/superBaseF2D.h"
#include "functors/superBaseF3D.h"



/// All OpenLB code is contained in this namespace.
namespace olb {

class BlockStructure2D;


/** A BlockLatticeReduction2D functor takes the data of a SupperLattice3D
 *  with targetDim=1 and provides an interface to BlockF2D. A 2D plane is
 *  given by a origin and two span vector u and v. Interpolation as well as
 *  the reduction of the plane to a cuboid (finite plane) is done in the
 *  constructor.
 *
 *  This class is not intended to be derived from.
 */
template< typename T, template <typename U> class DESCRIPTOR >
class BlockLatticeReduction2D final : public BlockDataF2D<T,T> {
private:
  /// Data fields which hold the reduced data
  BlockData2D<T,T>* _tmpBlockData;
  /// Functor which is reduced
  SuperLatticeF2D<T,DESCRIPTOR>& _f;
  /// Origin of the cuboid
  std::vector<T> _origin;
  /// number of voxel of the longest side,
  ///  default: _resolution=600,
  ///  off: _resolution=0 -> the smallest spacing of cuboidGeometry is chosen
  int _resolution;
  /// spacing
  T _h;
  /// Size of the cuboid
  int _nx;
  int _ny;
  /// Specific ostream for the classname in each line
  mutable OstreamManager clout;
public:
  /// Construction based on a functor f which is reduced to a 2D cuboid
  BlockLatticeReduction2D(SuperLatticeF2D<T,DESCRIPTOR>& f, int resolution=600);
  /// Destruction
  ~BlockLatticeReduction2D();
  /// Updates and writes the data to blockData
  void update();
  /// Attention: overload virtual function from class BlockLatticeF2D
  BlockStructure2D& getBlockStructure() override;
private:
  /// Updates _h, _nx, _ny such that the longest side is _resolution voxels long
  void updateToWantedResolution();
};


/** A BlockLatticeReduction3D functor takes the data of a SupperLattice3D
 *  with targetDim=1 and provides an interface to BlockF2D. A 2D plane is
 *  given by a origin and two span vector u and v. Interpolation as well as
 *  the reduction of the plane to a cuboid (finite plane) is done in the
 *  constructor.
 *
 *  This class is not intended to be derived from.
 */
template< typename T, template <typename U> class DESCRIPTOR >
class BlockLatticeReduction3D final : public BlockDataF2D<T,T> {
private:
  /// Data fields which hold the reduced data
  BlockData2D<T,T>* _tmpBlockData;
  /// Functor which is reduced
  SuperLatticeF3D<T,DESCRIPTOR>& _f;
  /// Origin of the cuboid
  Vector<T,3> _origin;
  /// number of voxel of the longest side,
  ///  default: _resolution=600,
  ///  off: _resolution=0 -> the smallest spacing of cuboidGeometry is chosen
  int _resolution;
  /// spacing
  T _h;
  /// Orthogonal vectors which span the cuboid, both are of lenght _h
  Vector<T,3> _u;
  Vector<T,3> _v;
  /// Size of the cuboid
  int _nx;
  int _ny;
  /// Specific ostream for the classname in each line
  mutable OstreamManager clout;
public:
  /// Construction based on a functor f which is reduced to a 2D cuboid (with origin and span vectors u and v)
  BlockLatticeReduction3D(SuperLatticeF3D<T,DESCRIPTOR>& f, Vector<T,3>& u,
                          Vector<T,3>& v, int resolution=600, T const origin[3]=nullptr);
  /// Construction based on a functor f which is reduced to a 2D cuboid (with origin and normal)
  //BlockLatticeReduction3D( SuperLatticeF3D<T,DESCRIPTOR>& f, T const normal[3], int resolution=600,
  //                         T const origin[3]=nullptr);
  /// Construction based on a functor f which is reduced to a 2D cuboid (with normal)
  BlockLatticeReduction3D(SuperLatticeF3D<T,DESCRIPTOR>& f, T const normalX,
                          T const normalY, T const normalZ, int resolution=600);
  /// Construction based on a functor f which is reduced to a 2D cuboid (with origin and normal)
  BlockLatticeReduction3D(SuperLatticeF3D<T,DESCRIPTOR>& f, T const normalX,
                          T const normalY, T const normalZ, T const originX,
                          T const originY, T const originZ, int resolution=600);
  /// Destruction
  ~BlockLatticeReduction3D();
  /// Access to the 2D grid blockData by cartCoord[0] = x, cartCoord[1] = y
  //  bool operator() (T output[], const int cartCoord[]) override;
  /// Updates and writes the data to blockData
  void update();
  /// Attention: overload virtual function from class BlockLatticeF2D
  BlockStructure2D& getBlockStructure() override;
private:
  /// Sets deafault as center of the cuboidGeometry or user definded origin
  void constructOrigin(T const origin[3]=NULL);
  /// Returns max possible distance
  int computeMaxLatticeDistance();
  /// Computes _origin and _nx, _ny such that the cuboid is right inside cuboidGeomety and not too big (sets _origin, _nx, _ny)
  void constructCuboid(int maxLatticeDistance);
  /// Updates _h, _nx, _ny, _nz, _u, _v such that the longest side is _resolution voxels long
  void updateToWantedResolution();
};

} // end namespace olb

#endif
