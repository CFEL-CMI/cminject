/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef SUPER_BASE_F_2D_H
#define SUPER_BASE_F_2D_H


#include "functors/genericF.h"
#include "functors/blockBaseF2D.h"
#include "communication/superStructure2D.h"
#include "core/superData2D.h"
#include "core/superLattice2D.h"
#include "core/units.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template<typename T, template<typename U> class Lattice> class SuperLattice2D;
template<typename T, typename BaseType> class SuperData2D;
template<typename T> class SuperStructure2D;
template<typename T> class BlockF2D;

/// represents all functors that operate on a SuperStructure in general
template <typename T, typename W = T>
class SuperF2D : public GenericF<W, int> {
protected:
  // ctor
  SuperF2D(SuperStructure2D<T>& superStructure, int targetDim);
  // dtor
  ~SuperF2D();
  // SuperFunctor consists of several BlockFuntors
  std::vector<BlockF2D<T>* > _blockF;
  SuperStructure2D<T>& _superStructure;
public:
  SuperF2D<T,W>& operator-(SuperF2D<T,W>& rhs);
  SuperF2D<T,W>& operator+(SuperF2D<T,W>& rhs);
  SuperF2D<T,W>& operator*(SuperF2D<T,W>& rhs);
  SuperF2D<T,W>& operator/(SuperF2D<T,W>& rhs);

  /// \return _superStructure
  SuperStructure2D<T>& getSuperStructure();
  /// \return _blockF[iCloc]
  BlockF2D<T>& getBlockF(int iCloc);
};


/// Functor from `SuperData2D`
template<typename T, typename BaseType>
class SuperDataF2D : public SuperF2D<T,BaseType> {
protected:
  /// `SuperData2D` object this functor was created from
  SuperData2D<T,BaseType>& _superData;
public:
  /// Constructor from `SuperData2D` - stores `_superData` reference
  SuperDataF2D(SuperData2D<T,BaseType>& superData);
  /// Operator for this functor - copies data from `_superData` object into output
  bool operator() (T output[], const int input[]);
  /// Getter for `_superData`
  SuperData2D<T,BaseType>& getSuperData();
};

/// identity functor for memory management
template <typename T, typename W>
class SuperIdentity2D final : public SuperF2D<T,W> {
protected:
  SuperF2D<T,W>& _f;
public:
  SuperIdentity2D(SuperF2D<T,W>& f);
  //  ~SuperLatticeIdentity2D();
  // access operator should not delete f, since f still has the identity as child
  bool operator() (W output[], const int input[]);
};

/// represents all functors that operate on a SuperLattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeF2D : public SuperF2D<T,T> {
protected:
  SuperLatticeF2D(SuperLattice2D<T,DESCRIPTOR>& superLattice, int targetDim);

  SuperLattice2D<T,DESCRIPTOR>& _sLattice;
public:
  SuperLattice2D<T,DESCRIPTOR>& getSuperLattice();
};

/// represents all functors that operate on a Lattice with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysF2D : public SuperLatticeF2D<T,DESCRIPTOR> {
protected:
  SuperLatticePhysF2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                      const LBconverter<T>& converter, int targetDim);
  const LBconverter<T>& _converter;
public:
  LBconverter<T> const& getConverter() const;
};


} // end namespace olb

#endif
