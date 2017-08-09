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

#ifndef SUPER_INDICATOR_F_3D_H
#define SUPER_INDICATOR_F_3D_H


#include "geometry/superGeometry3D.h"
#include "superIndicatorBaseF3D.h"


namespace olb {

template<typename T> class SuperGeometry3D;


/// Indicator Functor from Material Numbers
template <typename T>
class SuperIndicatorMaterial3D : public SuperIndicatorF3D<T> {
protected:
  SuperGeometry3D<T>& _superGeometry;
  std::vector<int> _materialNumbers;
public:
  SuperIndicatorMaterial3D (SuperGeometry3D<T>& rhs, std::vector<int> materialNumbers);
  virtual bool operator() (bool output[], const int input[]);
};


}

#endif
