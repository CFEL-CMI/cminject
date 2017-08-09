/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#include "superCalcF3D.h"
#include "superCalcF3D.hh"

namespace olb {

template class SuperCalc3D<double,double>;
template class SuperCalc3D<double,int>;
template class SuperCalc3D<double,bool>;

template class SuperPlus3D<double,double>;
template class SuperPlus3D<double,int>;
template class SuperPlus3D<double,bool>;

template class SuperMinus3D<double,double>;
template class SuperMinus3D<double,int>;
template class SuperMinus3D<double,bool>;

template class SuperMultiplication3D<double,double>;
template class SuperMultiplication3D<double,int>;
template class SuperMultiplication3D<double,bool>;

template class SuperDivision3D<double,double>;
template class SuperDivision3D<double,int>;

} // end namespace olb
