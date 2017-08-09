/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2015 Patrick Nathen, Mathias J. Krause
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

#ifndef TURBULENT_F_3D_H
#define TURBULENT_F_3D_H


#include "functors/blockBaseF3D.h"
#include "functors/superBaseF3D.h"
#include "functors/indicator/indicatorBaseF3D.h"


/** These are functors used for turbulent flows. Some like AMD have an execute member
 *  function which writes the data into the external field of a lattice descriptor.
 */

namespace olb {

/// functor to get pointwise yPlus from rho, shear stress and local density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeYplus3D : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  IndicatorF3D<T>&    _indicator;
  const int           _material;
public:
  SuperLatticeYplus3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter,
                      SuperGeometry3D<T>& superGeometry, IndicatorF3D<T>& indicator,
                      const int material );
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filtering on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
/*template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeADM3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  T _sigma;
  int _order;
  bool _adaptive;
  const LBconverter<T>& _converter;

public:
  BlockLatticeADM3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, T sigma, int order, bool adaptive, const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
  void execute(const int input[]);
  void execute();

  private:
  const  int _localAvDissBeginsAt = DESCRIPTOR<T>::ExternalField::localAvDissBeginsAt;
  const  int _localAvTKEBeginsAt = DESCRIPTOR<T>::ExternalField::localAvTKEBeginsAt;

};

/// functor to get pointwise ecplicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeADM3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  T _sigma;
  int _order;
  bool _adaptive;
  const LBconverter<T>& _converter;
public:
  SuperLatticeADM3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, T sigma, int order, bool adaptive, const LBconverter<T>& converter);
  bool operator() (T output[], const int input[]);
  void execute(SuperGeometry3D<T>& superGeometry, const int material);
};
*/
/// functor to get pointwise finite difference Dissipation on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysDissipationFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const LBconverter<T>& converter;
public:
  BlockLatticePhysDissipationFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& _converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDissipationFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  const LBconverter<T>& converter;
public:
  SuperLatticePhysDissipationFD3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& _converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise finite difference effective Dissipation on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysEffectiveDissipationFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const LBconverter<T>& converter;
  T smagoConst;
public:
  BlockLatticePhysEffectiveDissipationFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& _converter, T _smagoConst);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysEffectiveDissipationFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  const LBconverter<T>& converter;
  T smagoConst;
public:
  SuperLatticePhysEffectiveDissipationFD3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& _converter, T _smagoConst);
  bool operator() (T output[], const int input[]);
};

/*
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeSigmaADM3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:

  private:
  const  int _localSigmaADMBeginsAt = DESCRIPTOR<T>::ExternalField::localSigmaADMBeginsAt;
public:
  BlockLatticeSigmaADM3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeSigmaADM3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
public:
  SuperLatticeSigmaADM3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};
*/

} // end namespace olb

#endif


