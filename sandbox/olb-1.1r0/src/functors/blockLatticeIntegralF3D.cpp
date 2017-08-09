
/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#include "functors/blockLatticeIntegralF3D.h"
#include "functors/blockLatticeIntegralF3D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class BlockMax3D<double,descriptors::D3Q19Descriptor>;
template class BlockSum3D<double,descriptors::D3Q19Descriptor>;
template class BlockIntegral3D<double,descriptors::D3Q19Descriptor>;
template class BlockL1Norm3D<double,descriptors::D3Q19Descriptor>;
template class BlockL223D<double,descriptors::D3Q19Descriptor>;
template class BlockGeometryFaces3D<double>;
template class BlockLatticePhysDrag3D<double,descriptors::D3Q19Descriptor>;
template class BlockLatticePhysCorrDrag3D<double,descriptors::D3Q19Descriptor>;

template class BlockLatticeInterpPhysVelocity3D<double,descriptors::D3Q19Descriptor>;

} // end namespace olb


