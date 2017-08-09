/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ENTROPIC_LB_DYNAMICS_H
#define ENTROPIC_LB_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

template<typename T, template<typename U> class Lattice> class Cell;


/// Implementation of the entropic collision step
template<typename T, template<typename U> class Lattice>
class EntropicDynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  EntropicDynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual EntropicDynamics<T,Lattice>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
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
private:
  /// computes the entropy function H(f)=sum_i f_i*ln(f_i/t_i)
  T computeEntropy(const T f[]);
  /// computes the entropy growth H(f)-H(f-alpha*fNeq)
  T computeEntropyGrowth(const T f[], const T fNeq[], const T &alpha);
  /// computes the entropy growth derivative
  /// dH/dalpha=-sum_i fNeq_i*ln((f_i-alpha*fNeq_i)/t_i)
  T computeEntropyGrowthDerivative(const T f[], const T fNeq[], const T &alpha);
  /// Get the alpha parameter
  bool getAlpha(T &alpha, const T f[], const T fNeq[]);

  T omega;  ///< relaxation parameter
};

/// Implementation of the forced entropic collision step
template<typename T, template<typename U> class Lattice>
class ForcedEntropicDynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  ForcedEntropicDynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual ForcedEntropicDynamics<T,Lattice>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
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
private:
  /// computes the entropy function H(f)=sum_i f_i*ln(f_i/t_i)
  T computeEntropy(const T f[]);
  /// computes the entropy growth H(f)-H(f-alpha*fNeq)
  T computeEntropyGrowth(const T f[], const T fNeq[], const T &alpha);
  /// computes the entropy growth derivative
  /// dH/dalpha=-sum_i fNeq_i*ln((f_i-alpha*fNeq_i)/t_i)
  T computeEntropyGrowthDerivative(const T f[], const T fNeq[], const T &alpha);
  /// Get the alpha parameter
  bool getAlpha(T &alpha, const T f[], const T fNeq[]);

  T omega;  ///< relaxation parameter

  static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
  static const int sizeOfForce   = Lattice<T>::ExternalField::sizeOfForce;
};

}

#endif
