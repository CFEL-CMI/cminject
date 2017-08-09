/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#include "superCalcF2D.h"
#include "superCalcF2D.hh"

namespace olb {

template class SuperCalc2D<double,double>;
template class SuperCalc2D<double,int>;
template class SuperCalc2D<double,bool>;

template class SuperPlus2D<double,double>;
template class SuperPlus2D<double,int>;
template class SuperPlus2D<double,bool>;

template class SuperMinus2D<double,double>;
template class SuperMinus2D<double,int>;
template class SuperMinus2D<double,bool>;

template class SuperMultiplication2D<double,double>;
template class SuperMultiplication2D<double,int>;
template class SuperMultiplication2D<double,bool>;

template class SuperDivision2D<double,double>;
template class SuperDivision2D<double,int>;

} // end namespace olb
