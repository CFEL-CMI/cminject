/*
 *  Copyright (C) 2015 Marie-Luise Maier, Mathias J. Krause
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

/** Alberto Di Renzo, Francesco Paolo Di Maio:
 * "Comparison of contact-force models for the simulation of collisions in
 * DEM-based granular ow codes",
 * Chemical Engineering Science 59 (2004) 525 – 541
 */

#ifndef HERTZMINDLINDERESIEWICZ3D_H
#define HERTZMINDLINDERESIEWICZ3D_H

#include <cmath>
#include "functors/superLatticeLocalF3D.h"
#include "particles/particleSystem3D.h"
#include "forces.h"
#include "core/units.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, template<
           typename W> class DESCRIPTOR>
class HertzMindlinDeresiewicz3D: public Force3D<T, PARTICLETYPE> {

public:
  HertzMindlinDeresiewicz3D(T G1, T G2, T v1, T v2, T scale1 = T(1.), T scale2 = T(1.));
  virtual ~HertzMindlinDeresiewicz3D()
  {
  }
  ;
  void applyForce(typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
                  ParticleSystem3D<T, PARTICLETYPE> &pSys);

private:
  T _G1;  // Shear modulus (PA)
  T _G2;  // Shear modulus (PA)
  T _v1;  // poisson ratio
  T _v2;  // poisson ratio
  T _scale1; // scales normal value of force
  T _scale2; // scales tangential value of force
  T E1, E2; // E-Modul Particle
  T eE; // equivalent combined E-Modul
  T eG; // equivalent combined E-Modul
};

}

#endif
