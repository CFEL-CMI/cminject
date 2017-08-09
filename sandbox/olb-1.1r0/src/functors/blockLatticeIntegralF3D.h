/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_INTEGRAL_F_3D_H
#define BLOCK_LATTICE_INTEGRAL_F_3D_H


#include "functors/genericF.h"
#include "functors/blockBaseF3D.h"
#include "geometry/blockGeometry3D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


template<typename T, template<typename U> class Lattice> class BlockLattice3D;

/// BlockMax3D returns the max in each component of all points of a certain material
template <typename T, template <typename U> class DESCRIPTOR>
class BlockMax3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockMax3D(BlockLatticeF3D<T,DESCRIPTOR>& f, BlockGeometry3D<T>& blockGeometry,
             int material);
  bool operator() (T output[], const int input[]);
};


/// BlockSum3D sums over all cells of a certain material number
template <typename T, template <typename U> class DESCRIPTOR>
class BlockSum3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockSum3D(BlockLatticeF3D<T,DESCRIPTOR>& f, BlockGeometry3D<T>& blockGeometry,
             int material);
  bool operator() (T output[], const int input[]);
};


/// BlockIntegral3D
template <typename T, template <typename U> class DESCRIPTOR>
class BlockIntegral3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockIntegral3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                  BlockGeometry3D<T>& blockGeometry, int material);
  bool operator() (T output[], const int input[]);
};


/// BlockL1Norm3D returns componentwise the l1 norm
template <typename T, template <typename U> class DESCRIPTOR>
class BlockL1Norm3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockL1Norm3D(BlockLatticeF3D<T,DESCRIPTOR>& f, BlockGeometry3D<T>& blockGeometry, int material);
  bool operator() (T output[], const int input[]);
};


/// BlockL223D returns componentwise the squared l2-norm
template <typename T, template <typename U> class DESCRIPTOR>
class BlockL223D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockL223D(BlockLatticeF3D<T,DESCRIPTOR>& f,
             BlockGeometry3D<T>& blockGeometry,
             int material);
  bool operator() (T output[], const int input[]);
};


/// BlockGeometryFaces3D counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T>
class BlockGeometryFaces3D final : public GenericF<T,int> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  int _material;
  const LBconverter<T>& _converter;
public:
  BlockGeometryFaces3D(BlockGeometryStructure3D<T>& blockGeometry, int material,
                       const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/** functor to get pointwise phys force acting on a boundary with a given
 *  material on local lattice, if globIC is not on
 *  the local processor, the returned vector is empty
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysDrag3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysDrag3D(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                         BlockGeometry3D<T>& blockGeometry, int material,
                         const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/** functor to get pointwise phys force acting on a boundary with a given
 *  material on local lattice, if globIC is not on
 *  the local processor, the returned vector is empty
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysCorrDrag3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysCorrDrag3D(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                             BlockGeometry3D<T>& blockGeometry, int material,
                             const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};




} // end namespace olb

#endif
