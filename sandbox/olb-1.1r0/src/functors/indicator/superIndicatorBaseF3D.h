/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Benjamin Förster
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

#ifndef SUPER_INDICATOR_BASE_F_3D_H
#define SUPER_INDICATOR_BASE_F_3D_H

#include "functors/genericF.h"
#include "functors/superBaseF3D.h"
#include "communication/superStructure3D.h"


namespace olb {

template<typename T, typename W> class SuperF3D;
template<typename T> class SuperStructure3D;


/// Base indicator functor (discrete)
/**
 * Provides Union, Without and Intersection arithmetic.
 * _Note: `operator()` must be overloaded by child classes._
 */
template <typename T>
class SuperIndicatorF3D : public SuperF3D<T,bool> {
public:
  SuperIndicatorF3D(SuperStructure3D<T>& superStructure);
};


/// `SuperIndicatorF3D` from `IndicatorF3D`
template <typename T>
class SuperIndicatorFfromIndicatorF3D : public SuperIndicatorF3D<T> {
protected:
  IndicatorF3D<T>& _indicatorF;
public:
  SuperIndicatorFfromIndicatorF3D(IndicatorF3D<T>& indicatorF, SuperStructure3D<T>& superStructure);
  virtual bool operator() (bool output[], const int input[]);
};


}

#endif
