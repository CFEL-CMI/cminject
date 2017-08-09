/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2011 Jonas Latt, Mathias J. Krause,
 *  Jonas Kratzke
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

/** \file
 * Unit handling -- template instantiation.
 */

#include "units.h"
#include "units.hh"


namespace olb {

template class LBunits<double>;

template class LBconverter<double>;
//template class RTLBconverter<double>;
//template class RTLBConstconverter<double>;
//template class TLBconverter<double>;

template void writeLogFile
(LBunits<double> const& converter, std::string const& title);

template void writeLogFile
(LBconverter<double> const& converter, std::string const& title);

template LBconverter<double>* createLBconverter
(XMLreader const& params);

}
