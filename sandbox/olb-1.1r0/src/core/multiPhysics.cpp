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


#include "multiPhysics.h"
#include "dynamics/latticeDescriptors.h"
#include "dynamics/latticeDescriptors.hh"


namespace olb {

using namespace descriptors;

namespace multiPhysics {

template<>
MultiPhysicsId getMultiPhysicsScalarId<int>()
{
  return IntScalarFieldId;
}

template<>
MultiPhysicsId getMultiPhysicsScalarId<float>()
{
  return FloatScalarFieldId;
}

template<>
MultiPhysicsId getMultiPhysicsScalarId<double>()
{
  return DoubleScalarFieldId;
}


template<>
MultiPhysicsId getMultiPhysicsTensorId<int,2>()
{
  return IntTensorField2Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<float,2>()
{
  return FloatTensorField2Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<double,2>()
{
  return DoubleTensorField2Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<int,3>()
{
  return IntTensorField3Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<float,3>()
{
  return FloatTensorField3Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<double,3>()
{
  return DoubleTensorField3Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<int,4>()
{
  return IntTensorField4Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<float,4>()
{
  return FloatTensorField4Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<double,4>()
{
  return DoubleTensorField4Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<int,6>()
{
  return IntTensorField6Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<float,6>()
{
  return FloatTensorField6Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<double,6>()
{
  return DoubleTensorField6Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<int,9>()
{
  return IntTensorField9Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<float,9>()
{
  return FloatTensorField9Id;
}

template<>
MultiPhysicsId getMultiPhysicsTensorId<double,9>()
{
  return DoubleTensorField9Id;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<int,D2Q9Descriptor>()
{
  return IntD2Q9BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<float,D2Q9Descriptor>()
{
  return FloatD2Q9BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<double,D2Q9Descriptor>()
{
  return DoubleD2Q9BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<int,D3Q13Descriptor>()
{
  return IntD3Q13BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<float,D3Q13Descriptor>()
{
  return FloatD3Q13BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<double,D3Q13Descriptor>()
{
  return DoubleD3Q13BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<int,D3Q15Descriptor>()
{
  return IntD3Q15BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<float,D3Q15Descriptor>()
{
  return FloatD3Q15BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<double,D3Q15Descriptor>()
{
  return DoubleD3Q15BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<int,D3Q19Descriptor>()
{
  return IntD3Q19BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<float,D3Q19Descriptor>()
{
  return FloatD3Q19BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<double,D3Q19Descriptor>()
{
  return DoubleD3Q19BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<int,D3Q27Descriptor>()
{
  return IntD3Q27BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<float,D3Q27Descriptor>()
{
  return FloatD3Q27BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<double,D3Q27Descriptor>()
{
  return DoubleD3Q27BlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<int,ForcedD2Q9Descriptor>()
{
  return IntD2Q9WithForceBlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<float,ForcedD2Q9Descriptor>()
{
  return FloatD2Q9WithForceBlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<double,ForcedD2Q9Descriptor>()
{
  return DoubleD2Q9WithForceBlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<int,ForcedD3Q19Descriptor>()
{
  return IntD3Q19WithForceBlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<float,ForcedD3Q19Descriptor>()
{
  return FloatD3Q19WithForceBlockId;
}

template<>
MultiPhysicsId getMultiPhysicsBlockId<double,ForcedD3Q19Descriptor>()
{
  return DoubleD3Q19WithForceBlockId;
}

}

} // namespace olb
