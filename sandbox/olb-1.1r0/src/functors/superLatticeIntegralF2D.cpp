/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause,
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

#include "superLatticeIntegralF2D.h"
#include "superLatticeIntegralF2D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperMax2D<double>;
template class SuperMin2D<double>;
template class SuperSum2D<double>;
template class SuperSumIndicator2D<double>;
template class SuperIntegral2D<double>;
template class SuperL1Norm2D<double>;
template class SuperL2Norm2D<double>;
template class SuperLinfNorm2D<double>;
template class SuperL222D<double>;
template class SuperGeometryFaces2D<double>;
template class SuperGeometryFacesIndicator2D<double>;
template class SuperLatticePhysDrag2D<double,descriptors::D2Q9Descriptor>;
template class SuperLatticePhysDragIndicator2D<double,descriptors::D2Q9Descriptor>;
template class SuperLatticePhysCorrDrag2D<double,descriptors::D2Q9Descriptor>;
template class SuperLatticeFlux2D<double,descriptors::D2Q9Descriptor>;
template class SuperLatticePhysVelocityFlux2D<double,descriptors::D2Q9Descriptor>;
template class SuperLatticePhysPressureFlux2D<double,descriptors::D2Q9Descriptor>;

} // end namespace olb


