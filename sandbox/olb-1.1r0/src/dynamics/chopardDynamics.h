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

/** \file
 * BGK Dynamics with adjustable speed of sound -- header file.
 */
#ifndef CHOPARD_DYNAMICS_H
#define CHOPARD_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class ChopardDynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  ChopardDynamics(T vs2_, T omega_, Momenta<T,Lattice>& momenta_);
  ChopardDynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual ChopardDynamics<T,Lattice>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Initialize cell at equilibrium distribution
  virtual void iniEquilibrium(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d]);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local value of any parameter
  virtual T getParameter(int whichParameter) const;
  /// Set local value of any parameter
  virtual void setParameter(int whichParameter, T value);
  /// Set local speed of sound
  void setVs2(T vs2_);
  /// Get local speed of sound
  T    getVs2() const;
public:
  static T chopardBgkCollision (
    Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d], T vs2, T omega);
  static T chopardEquilibrium (
    int iPop, T rho, const T u[Lattice<T>::d], T uSqr, T vs2 );
private:
  T vs2;    ///< speed of sound
  T omega;  ///< relaxation parameter
};

}

#endif
