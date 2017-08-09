/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Mathias J. Krause, Cyril Masquelier, Benjamin Förster, Albert Mink
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

#include "indicCalcF2D.h"
#include "indicCalcF2D.hh"


namespace olb {

// arithmetic helper class for Indical 1d functors
template class IndicCalc1D<double>;
template class IndicPlus1D<double>;
template class IndicMinus1D<double>;
template class IndicMultiplication1D<double>;

// arithmetic helper class for Indical 2d functors
template class IndicCalc2D<double>;
template class IndicPlus2D<double>;
template class IndicMinus2D<double>;
template class IndicMultiplication2D<double>;

// arithmetic helper class for Indical 2d functors Smooth
template class SmoothIndicCalc2D<double,double>;
template class SmoothIndicPlus2D<double,double>;

}

