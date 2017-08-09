/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
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

#ifndef INTERPOLATION_F_3D_H
#define INTERPOLATION_F_3D_H


#include "functors/analyticalF.h"
#include "functors/blockBaseF3D.h"
#include "functors/superBaseF3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/blockGeometry3D.h"
#include "geometry/superGeometry3D.h"



namespace olb {

/// a class used to convert a block functor to analytical functor
template <typename T, typename W = T>
class AnalyticalFfromBlockF3D final : public AnalyticalF3D<T,W> {
protected:
  BlockF3D<W>& _f;
  Cuboid3D<T>& _cuboid;
public:
  AnalyticalFfromBlockF3D(BlockF3D<W>& f, Cuboid3D<T>& cuboid);
  bool operator() (W output[], const T physC[]);
};

/// a class used to convert a block functor to analytical functor
template <typename T, typename W = T>
class SpecialAnalyticalFfromBlockF3D final : public AnalyticalF3D<T,W> {
protected:
  BlockF3D<W>& _f;
  Cuboid3D<T>& _cuboid;
  Vector<T,3> _delta;
public:
  SpecialAnalyticalFfromBlockF3D(BlockF3D<W>& f, Cuboid3D<T>& cuboid, Vector<T,3> delta);
  bool operator() (W output[], const T physC[]);
};

/// a class used to convert super lattice functions to analytical functions
template <typename T, typename W = T>
class AnalyticalFfromSuperF3D final : public AnalyticalF3D<T,W> {
protected:
  SuperF3D<T,W>&                            _f;
  CuboidGeometry3D<T>&                      _cuboidGeometry;
  bool                                      _communicateToAll;
  int                                       _overlap;
  bool                                      _communicateOverlap;
  std::vector<AnalyticalFfromBlockF3D<T,W>* > _analyticalFfromBlockF;
public:
  AnalyticalFfromSuperF3D(SuperF3D<T,W>& f, bool communicateToAll=false, int overlap=-1,
                          bool communicateOverlap = true);
  virtual ~AnalyticalFfromSuperF3D();
  bool operator() (W output[], const T physC[]);
};


/**
 *  class used to convert analytical functions to lattice functions
 *  input functions are interpreted as SI->SI units, the resulting lattice
 *  function will map lattice->lattice units
 */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeFfromAnalyticalF3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  AnalyticalF3D<T,T>& _f;
public:
  SuperLatticeFfromAnalyticalF3D(AnalyticalF3D<T,T>& f,
                                 SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};


//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeFfromAnalyticalF3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  AnalyticalF3D<T,T>&  _f;
  BlockGeometry3D<T>&  _superGeometry;
  CuboidGeometry3D<T>& _cuboidGeometry;
public:
  BlockLatticeFfromAnalyticalF3D(AnalyticalF3D<T,T>& f,
                                 BlockLattice3D<T,DESCRIPTOR>& sLattice,
                                 BlockGeometry3D<T>& superGeometry,
                                 CuboidGeometry3D<T>& cuboidGeometry);
  bool operator() (T output[], const int input[]);
};

//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
class SmoothBlockIndicator3D final : public BlockDataF3D<T,T> {
protected:
  IndicatorF3D<T>&  _f;
  T _h;
public:
  SmoothBlockIndicator3D(IndicatorF3D<T>& f, T h);
  //bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeInterpPhysVelocity3Degree3D final : public
  BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  LBconverter<T>& _conv;
  Cuboid3D<T>* _cuboid;
  int _overlap;
  int _range;
public:
  BlockLatticeInterpPhysVelocity3Degree3D(
    BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
    LBconverter<T>& conv, Cuboid3D<T>* c, int overlap, int range);
  BlockLatticeInterpPhysVelocity3Degree3D(
    const BlockLatticeInterpPhysVelocity3Degree3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[], const int input[])
  {
    return false;
  }
  void operator() (T output[], const T input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeInterpPhysVelocity3Degree3D final : public
  SuperLatticeF3D<T,DESCRIPTOR> {
private:
  std::vector<BlockLatticeInterpPhysVelocity3Degree3D<T,DESCRIPTOR>* >
  bLattices;
public:
  SuperLatticeInterpPhysVelocity3Degree3D(
    SuperLattice3D<T,DESCRIPTOR>& sLattice, LBconverter<T>& conv,
    int range=1);
  bool operator() (T output[], const int input[])
  {
    return 0;
  }
  void operator()(T output[], const T input[], const int iC);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeInterpDensity3Degree3D final : public
  BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  LBconverter<T>& _conv;
  Cuboid3D<T>* _cuboid;
  int _overlap;
  int _range;
public:
  BlockLatticeInterpDensity3Degree3D(
    BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
    LBconverter<T>& conv, Cuboid3D<T>* c, int overlap, int range);
  BlockLatticeInterpDensity3Degree3D(
    const BlockLatticeInterpDensity3Degree3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[], const int input[])
  {
    return false;
  }
  void operator() (T output[], const T input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeInterpDensity3Degree3D final : public
  SuperLatticeF3D<T,DESCRIPTOR> {
private:
  std::vector<BlockLatticeInterpDensity3Degree3D<T,DESCRIPTOR>* > bLattices;
public:
  SuperLatticeInterpDensity3Degree3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     LBconverter<T>& conv, int range=1);
  bool operator() (T output[], const int input[])
  {
    return 0;
  }
  void operator()(T output[], const T input[], const int iC);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeSmoothDiracDelta3D final : public
  BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  Cuboid3D<T>* _cuboid;
  LBconverter<T>& _conv;
  BlockGeometryStructure3D<T>& _bGeometry;
  int _overlap;
public:
  BlockLatticeSmoothDiracDelta3D(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                 BlockGeometryStructure3D<T>& blockGeometry,
                                 Cuboid3D<T>* c, LBconverter<T>& conv, int overlap);
  BlockLatticeSmoothDiracDelta3D(
    const BlockLatticeSmoothDiracDelta3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[], const int input[])
  {
    return false;
  }
  void operator() (T delta[4][4][4], const T physPosP[3]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeSmoothDiracDelta3D final : public
  SuperLatticeF3D<T,DESCRIPTOR> {
private:
  std::vector<BlockLatticeSmoothDiracDelta3D<T,DESCRIPTOR>* > bLattices;
public:
  SuperLatticeSmoothDiracDelta3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 SuperGeometry3D<T>& superGeometry,
                                 LBconverter<T>& conv);
  void operator()(T delta[4][4][4], const T physPos[3], const int iC);
  bool operator()(T output[], const int input[])
  {
    return false;
  };
};



} // end namespace olb

#endif
