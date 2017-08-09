/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2015 Mathias J. Krause, Jonas Latt, Patrick Nathen
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
 * BGK Dynamics with adjusted omega -- header file.
 */
#ifndef SMAGORINSKY_BGK_DYNAMICS_H
#define SMAGORINSKY_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "core/units.h"
#include "core/units.hh"

namespace olb {


/// Implementation of the ADM BGK collision step

template<typename T, template<typename U> class Lattice>
class ADMBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  ADMBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
private:
  T omega;
  static const int rhoIsAt = Lattice<T>::ExternalField::rhoIsAt;
  static const int velocityBeginsAt = Lattice<T>::ExternalField::velocityBeginsAt;
  static const int sizeOfVelocity = Lattice<T>::ExternalField::sizeOfVelocity;
};

/// Implementation of the ForcedADMBGK collision step
template<typename T, template<typename U> class Lattice>
class ForcedADMBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  ForcedADMBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_);

  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
private:
  T omega;
  static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
  static const int sizeOfForce   = Lattice<T>::ExternalField::sizeOfForce;
  static const int filRhoIsAt = Lattice<T>::ExternalField::filRhoIsAt;
  static const int localFilVelXBeginsAt = Lattice<T>::ExternalField::localFilVelXBeginsAt;
  static const int localFilVelYBeginsAt = Lattice<T>::ExternalField::localFilVelYBeginsAt;
  static const int localFilVelZBeginsAt = Lattice<T>::ExternalField::localFilVelZBeginsAt;
};


/// Implementation of the consistent Smagorinsky BGK collision step
///
/// Consistent subgrid scale modelling for lattice Boltzmann methods
/// Orestis Malaspinas and Pierre Sagaut
/// Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
/// DOI: http://dx.doi.org/10.1017/jfm.2012.155

template<typename T, template<typename U> class Lattice>
class ConSmagorinskyBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  ConSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_,
                            T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// should be remove --> David
  T computeOmega(T omega0, T preFactor_, T rho, T pi[util::TensorVal<Lattice<T> >::n] );

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};


/// Implementation of the consistent Strain Smagorinsky BGK collision step
///
/// Consistent subgrid scale modelling for lattice Boltzmann methods
/// Orestis Malaspinas and Pierre Sagaut
/// Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
/// DOI: http://dx.doi.org/10.1017/jfm.2012.155


template<typename T, template<typename U> class Lattice>
class ConStrainSmagorinskyBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  ConStrainSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
                                  T smagoConst_=T(.1), T dx_=T(1), T dt_=T(1));
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFactor_, T rho, T pi[util::TensorVal<Lattice<T> >::n]);
  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};


/// Implementation of a the dynamic Smarorinsky BGK collision step
template<typename T, template<typename U> class Lattice>
class DynSmagorinskyBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  DynSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFactor_, T rho, T pi[util::TensorVal<Lattice<T> >::n],
                 Cell<T,Lattice>& cell);
  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  T dx;
  T dt;
  const static T smagoConstIsAt = Lattice<T>::ExternalField::smagoConstIsAt;
};


/// Implementation of a Shear Smarorinsky BGK collision step
/// Shown good results for wall-bounded flows
/// Leveque et al.: Shear-Improved Smagorinsky Model for Large-Eddy Simulation
/// of Wall-Bounded Turbulent Flows
/// DOI: http://dx.doi.org/10.1017/S0022112006003429

template<typename T, template<typename U> class Lattice>
class ShearSmagorinskyBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  ShearSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_,
                              T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_, int iT);

private:

  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0_, T rho_, T pi_[util::TensorVal<Lattice<T> >::n],
                 Cell<T,Lattice>& cell, int iT);

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Pre factor speeds up the collision step
  T preFactor;
  T dx;
  T dt;
  const static int avShearIsAt = Lattice<T>::ExternalField::avShearIsAt;
};

/// Implementation of the ForcedBGK collision step
template<typename T, template<typename U> class Lattice>
class ShearSmagorinskyForcedBGKdynamics : public ForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ShearSmagorinskyForcedBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
                                    T smagoConst_, T dx_ = 1, T dt_ = 1);
  ///collids
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);

  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_, int iT);

private:

  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0_, T rho_, T pi_[util::TensorVal<Lattice<T> >::n],
                 Cell<T,Lattice>& cell, int iT);

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Pre factor speeds up the collision step
  T preFactor;
  T dx;
  T dt;
  const static int avShearIsAt = Lattice<T>::ExternalField::avShearIsAt;
};


/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class SmagorinskyBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_,
                         T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFactor_, T rho,
                 T pi[util::TensorVal<Lattice<T> >::n] );

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};


/// Implementation of the ForcedBGK collision step
template<typename T, template<typename U> class Lattice>
class SmagorinskyForcedBGKdynamics : public ForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyForcedBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
                               T smagoConst_, T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFactor_, T rho,
                 T pi[util::TensorVal<Lattice<T> >::n]);

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};


/// Implementation of the ForcedBGK collision step
template<typename T, template<typename U> class Lattice>
class SmagorinskyLinearVelocityForcedBGKdynamics : public ForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyLinearVelocityForcedBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
      T smagoConst_, T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFactor_, T rho,
                 T pi[util::TensorVal<Lattice<T> >::n]);

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};

/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class KrauseBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  KrauseBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_,
                    T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  void computeOmega(T omega0, Cell<T,Lattice>& cell, T preFactor_, T rho,
                    T u[Lattice<T>::d],
                    T newOmega[Lattice<T>::q] );

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};


}

#endif
