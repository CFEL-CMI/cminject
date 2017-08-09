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

#ifndef BLOCK_LATTICE_INTEGRAL_F_2D_H
#define BLOCK_LATTICE_INTEGRAL_F_2D_H


#include "functors/genericF.h"
#include "functors/blockBaseF2D.h"
#include "geometry/blockGeometry2D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


template<typename T, template<typename U> class Lattice> class BlockLattice2D;


/// BlockMax2D returns the max in each component of all points of a certain material
template <typename T, template <typename U> class DESCRIPTOR>
class BlockMax2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockMax2D(BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry2D<T>& blockGeometry,
             int material);
  bool operator() (T output[], const int input[]);
};


/// BlockSum2D sums over all cells of a certain material number
template <typename T, template <typename U> class DESCRIPTOR>
class BlockSum2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockSum2D(BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry2D<T>& blockGeometry,
             int material);
  bool operator() (T output[], const int input[]);
};


/// BlockIntegral2D
template <typename T, template <typename U> class DESCRIPTOR>
class BlockIntegral2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockIntegral2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                  BlockGeometry2D<T>& blockGeometry, int material);
  bool operator() (T output[], const int input[]);
};


/// BlockL1Norm2D returns componentwise the l1 norm
template <typename T, template <typename U> class DESCRIPTOR>
class BlockL1Norm2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockL1Norm2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                BlockGeometry2D<T>& blockGeometry, int material);
  bool operator() (T output[], const int input[]);
};


/// BlockL222D returns componentwise the squared l2-norm
template <typename T, template <typename U> class DESCRIPTOR>
class BlockL222D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockL222D(BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry2D<T>& blockGeometry,
             int material);
  bool operator() (T output[], const int input[]);
};


/// functor counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T>
class BlockGeometryFaces2D final : public GenericF<T,int> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  int _material;
  const LBconverter<T>& _converter;
public:
  BlockGeometryFaces2D(BlockGeometryStructure2D<T>& blockGeometry, int material,
                       const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor counts to get the discrete surface for a smooth indicator circle in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)
/// and total surface, then it converts it into phys units
template <typename T>
class BlockGeometryFacesIndicator2D final : public GenericF<T,int> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  SmoothIndicatorCircle2D<T,T>& _indicator;
  int _material;
  const LBconverter<T>& _converter;
public:
  BlockGeometryFacesIndicator2D(BlockGeometryStructure2D<T>& blockGeometry, SmoothIndicatorCircle2D<T,T>& indicator,
                                int material, const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};



/** functor to get pointwise phys force acting on a boundary with a given
 *  material on local lattice, if globIC is not on
 *  the local processor, the returned vector is empty
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysDrag2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysDrag2D(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                         BlockGeometry2D<T>& blockGeometry, int material,
                         const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};


/** functor to get pointwise phys force acting on a boundary with a given
 *  material on local lattice, if globIC is not on
 *  the local processor, the returned vector is empty
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysCorrDrag2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysCorrDrag2D(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                             BlockGeometry2D<T>& blockGeometry, int material,
                             const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
};




} // end namespace olb

#endif
