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

#ifndef HERTZMINDLINDERESIEWICZ3D_HH
#define HERTZMINDLINDERESIEWICZ3D_HH

#include <cmath>

#include "hertzMindlinDeresiewicz3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, template<
           typename W> class DESCRIPTOR>
HertzMindlinDeresiewicz3D<T, PARTICLETYPE, DESCRIPTOR>::HertzMindlinDeresiewicz3D(
  T G1, T G2, T v1, T v2, T scale1, T scale2) :
  Force3D<T, PARTICLETYPE>(), _G1(G1), _G2(G2), _v1(v1), _v2(v2), _scale1(
    scale1), _scale2(scale2)
{
  // E-Modul Particle
  E1 = 2 * (1 + _v1) * _G1;
  E2 = 2 * (1 + _v2) * _G2;

  // equivalent combined E-Modul
  eE = (1 - pow(_v1, 2)) / E1 + (1 - pow(_v2, 2)) / E2;
  eE = 1 / eE;

  // equivalent combined E-Modul
  eG = (2.0 - _v1) / _G1 + (2 - _v2) / _G2;
  eG = 1. / eG;
}

template<typename T, template<typename U> class PARTICLETYPE, template<
           typename W> class DESCRIPTOR>
void HertzMindlinDeresiewicz3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& pSys)
{

  std::vector<std::pair<size_t, T>> ret_matches;
  // kind of contactDetection has to be chosen in application
  pSys.getContactDetection()->getMatches(pInt, ret_matches);

  const PARTICLETYPE<T>* p2 = NULL;

  typename std::vector<std::pair<size_t, T> >::iterator it =
    ret_matches.begin();

  // iterator walks through number of neighbored particles = ret_matches
  for (; it != ret_matches.end(); it++) {
    if (it->second != 0) {
      p2 = &pSys[it->first];
      // overlap
      T delta = (p2->getRad() + p->getRad()) - sqrt(it->second);

      // equivalent mass
      T M = p->getMass() * p2->getMass() / (p->getMass() + p2->getMass());
      // equivalent radius
      T R = p->getRad() * p2->getRad() / (p->getRad() + p2->getRad());
      // relative velocity
      std::vector < T > _velR(3, T());
      _velR[0] = -(p2->getVel()[0] - p->getVel()[0]);
      _velR[1] = -(p2->getVel()[1] - p->getVel()[1]);
      _velR[2] = -(p2->getVel()[2] - p->getVel()[2]);

      std::vector < T > _d(3, T());
      std::vector < T > _normal(3, T());

      //_d: vector from particle1 to particle2
      _d[0] = p2->getPos()[0] - p->getPos()[0];
      _d[1] = p2->getPos()[1] - p->getPos()[1];
      _d[2] = p2->getPos()[2] - p->getPos()[2];

      if (util::norm(_d) != T(0)) {
        _normal = util::normalize(_d);
      } else {
        return;
      }

      Vector<T, 3> d_(_d);
      Vector<T, 3> velR_(_velR);
      T dot = velR_[0]*_normal[0] + velR_[1]*_normal[1]+ velR_[2]*_normal[2];

      // normal part of relative velocity
      // normal relative to surface of particles at contact point
      std::vector < T > _velN(3, T());
      _velN[0] = dot * _normal[0];
      _velN[1] = dot * _normal[1];
      _velN[2] = dot * _normal[2];

      // tangential part of relative velocity
      // tangential relative to surface of particles at contact point
      std::vector < T > _velT(3, T());
      _velT[0] = _velR[0] - _velN[0];
      _velT[1] = _velR[1] - _velN[1];
      _velT[2] = _velR[2] - _velN[2];

      if (pow(p2->getRad() + p->getRad(), 2) > it->second) { // overlap

        // Force normal
        // spring constant in normal direction
        // (Alberto Di Renzo∗ , Francesco Paolo Di Maio, Chemical Engineering Science 59 (2004) 525 – 541)
        // Konstante kn aus (H. Kruggel-Endem, Powder Technology 171 (2007) 157-173) Table 3 / Brass Partikel
        T kn = 4 / 3. * sqrt(R) * eE;

        // part of mechanical force of spring in normal direction
        // Hertz Contact (P. A. Langston, Powder Technology 85 (1995))
        std::vector < T > Fs_n(3, T());
        Fs_n[0] = -kn * pow(delta, 1.5) * _normal[0];
        Fs_n[1] = -kn * pow(delta, 1.5) * _normal[1];
        Fs_n[2] = -kn * pow(delta, 1.5) * _normal[2];

        // part of mechanical force of damper in normal direction
        // damped linear spring (Cundall, Strack 1979)
        // (K.W. Chu, A.B. Yu, Powder Technology 179 (2008) 104 – 114)
        // damper constant in normal direction
        // Konstante eta_n aus (H. Kruggel-Endem, Powder Technology 171 (2007) 157-173) Table 3 / Brass Partikel
        T eta_n = 0.3 * sqrt(4.5 * M * sqrt(delta) * kn);

        std::vector < T > Fd_n(3, T());
        Fd_n[0] = -eta_n * _velN[0]*sqrt(delta);
        Fd_n[1] = -eta_n * _velN[1]*sqrt(delta);
        Fd_n[2] = -eta_n * _velN[2]*sqrt(delta);

        std::vector < T > F_n(3, T());
        F_n[0] = Fs_n[0] + Fd_n[0];
        F_n[1] = Fs_n[1] + Fd_n[1];
        F_n[2] = Fs_n[2] + Fd_n[2];

        // Force tangential
        // spring constant in tangential direction
        // (N.G. Deen, Chemical Engineering Science 62 (2007) 28 – 44)
        T kt = 2 * sqrt(2 * R) * _G1 / (2 - _v1) * pow(delta, 0.5);

        // damper constant in normal direction
        T eta_t = 2 * sqrt(2. / 7. * M * kt);

        // part of mechanical force of damper in tangential direction
        std::vector < T > F_t(3, T());
        F_t[0] = -eta_t * _velT[0];
        F_t[1] = -eta_t * _velT[1];
        F_t[2] = -eta_t * _velT[2];

        // entire force
        // factor _scale to prevent instability
        p->getForce()[0] += _scale1 * F_n[0] + _scale2 * F_t[0];
        p->getForce()[1] += _scale1 * F_n[1] + _scale2 * F_t[1];
        p->getForce()[2] += _scale1 * F_n[2] + _scale2 * F_t[2];
      }
    }
  }
}

}

#endif
