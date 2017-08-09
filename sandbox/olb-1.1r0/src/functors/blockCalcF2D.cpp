/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Albert Mink, Lukas Baron, Mathias J. Krause
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

#include "functors/blockCalcF2D.h"
#include "functors/blockCalcF2D.hh"

namespace olb {

template class BlockCalc2D<int>;
template class BlockCalc2D<double>;

template class BlockPlus2D<int>;
template class BlockPlus2D<double>;

template class BlockMinus2D<int>;
template class BlockMinus2D<double>;

template class BlockMultiplication2D<int>;
template class BlockMultiplication2D<double>;

template class BlockDivision2D<int>;
template class BlockDivision2D<double>;

} // end namespace olb
