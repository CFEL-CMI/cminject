/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2016 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Benjamin Förster
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

#include "superLatticeIntegralF3D.h"
#include "superLatticeIntegralF3D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperMin3D<double,double>;
template class SuperMin3D<double,int>;
template class SuperMin3D<double,bool>;

template class SuperMax3D<double,double>;
template class SuperMax3D<double,int>;
template class SuperMax3D<double,bool>;

template class SuperSum3D<double,double>;
template class SuperSum3D<double,int>;

template class SuperSumIndicator3D<double,double>;

template class SuperAverage3D<double,double>;

template class SuperIntegral3D<double,double>;
template class SuperIntegral3D<double,int>;

template class SuperLpNorm3D<double,double>;
template class SuperLpNorm3D<double,int>;
template class SuperL1Norm3D<double,double>;
template class SuperL1Norm3D<double,int>;
template class SuperL2Norm3D<double,double>;
template class SuperL2Norm3D<double,int>;
template class SuperLinfNorm3D<double,double>;
template class SuperLinfNorm3D<double,int>;

template class SuperGeometryFaces3D<double>;

template class SuperLatticePhysDrag3D<double,descriptors::D3Q19Descriptor>;
template class SuperLatticePhysDragIndicator3D<double,descriptors::D3Q19Descriptor>;
template class SuperLatticePhysCorrDrag3D<double,descriptors::D3Q19Descriptor>;
template class SuperLatticeFlux3D<double,descriptors::D3Q19Descriptor>;
template class SuperLatticePhysPressureFlux3D<double,descriptors::D3Q19Descriptor>;
template class SuperLatticePhysVelocityFlux3D<double,descriptors::D3Q19Descriptor>;

} // end namespace olb


