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

#ifndef SUPER_INDICATOR_F_3D_HH
#define SUPER_INDICATOR_F_3D_HH

#include "superIndicatorF3D.h"
#include "core/util.h"
#include<numeric>


namespace olb {


template <typename S>
SuperIndicatorMaterial3D<S>::SuperIndicatorMaterial3D (SuperGeometry3D<S>& rhs, std::vector<int> materialNumbers)
  : SuperIndicatorF3D<S>(rhs), _superGeometry(rhs), _materialNumbers(materialNumbers)
{
  std::string matString = std::accumulate(_materialNumbers.begin()+1, _materialNumbers.end(), std::to_string(_materialNumbers[0]),
  [](const std::string& a, int b) {
    return a + '_' + std::to_string(b);
  });
  this->getName() = "SuperIndicator_on_Material_" + matString;
}

template <typename S>
bool SuperIndicatorMaterial3D<S>::operator() (bool output[], const int input[])
{
  // iterate over material numbers and check if given point has that material number
  output[0] = false;
  for (int& m : _materialNumbers) {
    output[0] |= ( _superGeometry.get(input[0], input[1], input[2], input[3]) == m );
  }
  return output[0];
}


} // namespace olb

#endif
