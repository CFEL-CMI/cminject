/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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

#ifndef ZOU_HE_BOUNDARY_3D_H
#define ZOU_HE_BOUNDARY_3D_H

#include "boundaryCondition3D.h"

namespace olb {

////////// Factory function for Zou/He BC ///////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,Lattice>*
createZouHeBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block);

template<typename T, template<typename U> class Lattice>
OnLatticeBoundaryCondition3D<T,Lattice>*
createZouHeBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return createZouHeBoundaryCondition3D<T,Lattice,BGKdynamics<T,Lattice> >(block);
}



}


#endif
