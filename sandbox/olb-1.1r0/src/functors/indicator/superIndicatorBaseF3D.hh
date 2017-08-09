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

#ifndef SUPER_INDICATOR_BASE_F_3D_HH
#define SUPER_INDICATOR_BASE_F_3D_HH


#include "superIndicatorBaseF3D.h"

namespace olb {

template <typename T>
SuperIndicatorF3D<T>::SuperIndicatorF3D(SuperStructure3D<T>& superStructure)
  : SuperF3D<T, bool>(superStructure, 1)
{ }


template <typename T>
SuperIndicatorFfromIndicatorF3D<T>::SuperIndicatorFfromIndicatorF3D(IndicatorF3D<T>& indicatorF,
    SuperStructure3D<T>& superStructure)
  : SuperIndicatorF3D<T> (superStructure),
    _indicatorF(indicatorF)
{
  this->getName() = "SuperIndicator_from_" + _indicatorF.getName();
}


template <typename T>
bool SuperIndicatorFfromIndicatorF3D<T>::operator() (bool output[], const int input[])
{
  T physR[3];
  this->_superStructure.getCuboidGeometry().getPhysR(physR, input);
  _indicatorF(output, physR);
  return true;
}


} // namespace olb

#endif
