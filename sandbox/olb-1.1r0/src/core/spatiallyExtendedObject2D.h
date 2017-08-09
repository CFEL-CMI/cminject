/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Jonas Latt, Orestis Malaspinas
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


#ifndef SPATIALLY_EXTENDED_OBJECT_2D_H
#define SPATIALLY_EXTENDED_OBJECT_2D_H

#include "multiPhysics.h"

namespace olb {

struct SpatiallyExtendedObject2D {
  virtual ~SpatiallyExtendedObject2D() { }
  virtual SpatiallyExtendedObject2D* getComponent(int iBlock) =0;
  virtual SpatiallyExtendedObject2D const* getComponent(int iBlock) const =0;
  virtual multiPhysics::MultiPhysicsId getMultiPhysicsId() const =0;
};

} // namespace olb

#endif
