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


#ifndef MULTI_PHYSICS_H
#define MULTI_PHYSICS_H

namespace olb {

namespace multiPhysics {

typedef enum {
  UndefinedId           = 0,

  IntScalarFieldId      = 10000,
  FloatScalarFieldId    = 10001,
  DoubleScalarFieldId   = 10002,

  IntTensorField2Id     = 20020,
  FloatTensorField2Id   = 20021,
  DoubleTensorField2Id  = 20022,

  IntTensorField3Id     = 20030,
  FloatTensorField3Id   = 20031,
  DoubleTensorField3Id  = 20032,

  IntTensorField4Id     = 20040,
  FloatTensorField4Id   = 20041,
  DoubleTensorField4Id  = 20042,

  IntTensorField6Id     = 20060,
  FloatTensorField6Id   = 20061,
  DoubleTensorField6Id  = 20062,

  IntTensorField9Id     = 20090,
  FloatTensorField9Id   = 20091,
  DoubleTensorField9Id  = 20092,

  IntD2Q5BlockId        = 32040,
  FloatD2Q5BlockId      = 32041,
  DoubleD2Q5BlockId     = 32042,

  IntD2Q9BlockId        = 32090,
  FloatD2Q9BlockId      = 32091,
  DoubleD2Q9BlockId     = 32092,

  IntD3Q7BlockId        = 33070,
  FloatD3Q7BlockId      = 33071,
  DoubleD3Q7BlockId     = 33072,

  IntD3Q13BlockId       = 33130,
  FloatD3Q13BlockId     = 33131,
  DoubleD3Q13BlockId    = 33132,

  IntD3Q15BlockId       = 33150,
  FloatD3Q15BlockId     = 33151,
  DoubleD3Q15BlockId    = 33152,

  IntD3Q19BlockId       = 33190,
  FloatD3Q19BlockId     = 33191,
  DoubleD3Q19BlockId    = 33192,

  IntD3Q27BlockId       = 33270,
  FloatD3Q27BlockId     = 33271,
  DoubleD3Q27BlockId    = 33272,

  IntD2Q5WithForceBlockId    = 42040,
  FloatD2Q5WithForceBlockId  = 42041,
  DoubleD2Q5WithForceBlockId = 42042,

  IntD2Q9WithForceBlockId    = 42090,
  FloatD2Q9WithForceBlockId  = 42091,
  DoubleD2Q9WithForceBlockId = 42092,

  IntD3Q7WithForceBlockId    = 43070,
  FloatD3Q7WithForceBlockId  = 43071,
  DoubleD3Q7WithForceBlockId = 43072,

  IntD3Q13WithForceBlockId    = 43130,
  FloatD3Q13WithForceBlockId  = 43131,
  DoubleD3Q13WithForceBlockId = 43132,

  IntD3Q15WithForceBlockId    = 43150,
  FloatD3Q15WithForceBlockId  = 43151,
  DoubleD3Q15WithForceBlockId = 43152,

  IntD3Q19WithForceBlockId    = 43190,
  FloatD3Q19WithForceBlockId  = 43191,
  DoubleD3Q19WithForceBlockId = 43192,

  IntD3Q27WithForceBlockId    = 43270,
  FloatD3Q27WithForceBlockId  = 43271,
  DoubleD3Q27WithForceBlockId = 43272
}
MultiPhysicsId;

template<typename T>
MultiPhysicsId getMultiPhysicsScalarId()
{
  return UndefinedId;
}

template<typename T, int n>
MultiPhysicsId getMultiPhysicsTensorId()
{
  return UndefinedId;
}

template<typename T, template<typename U> class Lattice>
MultiPhysicsId getMultiPhysicsBlockId()
{
  return UndefinedId;
}

} // namespace multiPhysics

} // namespace olb

#endif
