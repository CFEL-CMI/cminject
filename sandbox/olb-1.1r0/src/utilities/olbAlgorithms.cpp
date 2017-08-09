/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Jonas Latt
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

#include "olbAlgorithms.h"

namespace olb {

namespace algorithm {

std::vector<int> primeFactor(int value)
{
  std::vector<int> primeFactors;
  int testFactor = 2;
  while (testFactor <= value) {
    if (value%testFactor==0) {
      value /= testFactor;
      primeFactors.push_back(testFactor);
    } else {
      ++testFactor;
    }
  }
  return primeFactors;
}

std::vector<int> evenRepartition(int value, int d)
{
  std::vector<int> primeFactors = primeFactor(value);
  std::vector<int> repartition(d);
  for (int iRep=0; iRep<d; ++iRep) {
    repartition[iRep] = 1;
  }
  int iDim=0;
  for (int iPrime=(int)(primeFactors.size()-1); iPrime>=0; --iPrime) {
    repartition[iDim] *= primeFactors[iPrime];
    iDim = (iDim+1)%d;
  }
  return repartition;
}

}

}
