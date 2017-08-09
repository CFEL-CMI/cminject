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

#ifndef INTERPOLATION_F_2D_H
#define INTERPOLATION_F_2D_H


#include "functors/analyticalF.h"
#include "functors/superBaseF2D.h"
#include "geometry/superGeometry2D.h"
#include "geometry/cuboidGeometry2D.h"

namespace olb {


template< typename T, template <typename U> class DESCRIPTOR> class SuperLattice2D;


/// converts lattice functions to analytical functions
template <typename T, template <typename U> class DESCRIPTOR>
class AnalyticalFfromSuperLatticeF2D final : public AnalyticalF2D<T,T> {
protected:
  SuperLatticeF2D<T,DESCRIPTOR>&  _f;
  CuboidGeometry2D<T>&            _cg;
  bool                            _communicateToAll;
  int                             _overlap;
public:
  AnalyticalFfromSuperLatticeF2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
                                 bool communicateToAll = false, int overlap = -1);
  bool operator() (T output[], const T physC[]);
};


/**
 *  a class used to convert analytical functions to lattice functions
 *  input functions are interpreted as SI->SI units, the resulting lattice
 *  function will map lattice->lattice units
 */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeFfromAnalyticalF2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
protected:
  AnalyticalF2D<T,T>&  _f;
  SuperGeometry2D<T>&  _sg;
public:
  SuperLatticeFfromAnalyticalF2D(AnalyticalF2D<T,T>& f, SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                 SuperGeometry2D<T>& sg);
  bool operator() (T output[], const int input[]);
};




} // end namespace olb

#endif
