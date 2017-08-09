/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2016 Jonas Latt, Mathias J. Krause
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

#ifndef SERIALIZER_IO_H
#define SERIALIZER_IO_H

#include "core/serializer.h"
#include <ostream>
#include <fstream>

namespace olb {

class Serializer;

/// processes data from a serializer to a given ostr, always in parallel
void serializer2ostr(Serializer& serializer, std::ostream& ostr, bool enforceUint=false);
/// processes an istr to a serializer, always in parallel
void istr2serializer(Serializer& serializer, std::istream& istr, bool enforceUint=false);

} // namespace olb

#endif
