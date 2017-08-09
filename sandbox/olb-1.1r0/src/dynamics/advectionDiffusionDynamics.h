/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008,2015 Orestis Malaspinas, Andrea Parmigiani, Albert Mink
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
#ifndef ADVECTION_DIFFUSION_DYNAMICS_H
#define ADVECTION_DIFFUSION_DYNAMICS_H

#include "advectionDiffusionLatticeDescriptors.h"
#include "dynamics/dynamics.h"

namespace olb {

// ========= the RLB advection diffusion dynamics ========//
/// it uses the regularized approximation that can be found in the thesis of J. Latt (2007).
template<typename T, template<typename U> class Lattice>
class AdvectionDiffusionRLBdynamics : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  AdvectionDiffusionRLBdynamics( T omega_, Momenta<T, Lattice>& momenta_ );
  /// Clone the object on its dynamic type.
  virtual AdvectionDiffusionRLBdynamics<T, Lattice>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const;
  /// Collision step
  virtual void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics );
  /// Collide with fixed velocity
  virtual void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
                              LatticeStatistics<T>& statistics );
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega( T omega_ );
private:
  T omega;  ///< relaxation parameter
};

// ========= the BGK advection diffusion dynamics ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, template<typename U> class Lattice>
class AdvectionDiffusionBGKdynamics : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  AdvectionDiffusionBGKdynamics( T omega, Momenta<T, Lattice>& momenta );
  /// Clone the object on its dynamic type.
  virtual AdvectionDiffusionBGKdynamics<T, Lattice>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const;
  /// Collision step
  virtual void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics );
  /// Collide with fixed velocity
  virtual void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
                              LatticeStatistics<T>& statistics );
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega( T omega );
private:
  T _omega;  ///< relaxation parameter
};

/**
 * Solves diffusion equation with additional sink term, according to Albert Mink et al 2016.
 *
 * \f[ \Delta \Phi = \frac{\sigma_a}{D} * \Phi \f]
 * where diffusion coefficient D is given by:
 * \f$ D = \frac{1}{3(\sigma_a + \sigma_s)} \f$
 * absorption and scattering coefficient:
 * \f$ \sigma_a \f$ and \f$ \sigma_s \f$
 *
 * \param omega       is relaxation parameter and always set to 1, since the above equation is divided by diffusionsCoeff.
 * \param sink        corresponds to the sink factor (grid independent), given by \f$ \frac{\sigma_a}{8*D} * N^{-2} \f$
 */
template<typename T, template<typename U> class Lattice>
class AdDiSinkBGKdynamics : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  AdDiSinkBGKdynamics( T omega, Momenta<T, Lattice>& momenta, T sink );
  /// Clone the object on its dynamic type.
  virtual AdDiSinkBGKdynamics<T, Lattice>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const;
  /// Collision step
  virtual void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics );
  /// Collide with fixed velocity
  virtual void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
                              LatticeStatistics<T>& statistics );
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega( T omega );
private:
  T _omega;
  T _sink;
};


/**
 * Solves RTE according Christopher McHardy et al 2016.
 * absorption and scattering coefficient:
 * \f$ \sigma_a \f$ and \f$ \sigma_s \f$
 *
 * \param omega             change into beta the extinction coefficient
 * \param singleScatAlbedo  is the single scattering albedo, given by \f$ \frac{\sigma_s}{sigma_a + sigma_s} \f$
 */
template<typename T, template<typename U> class Lattice>
class AdDiAnisoBGKdynamics : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  AdDiAnisoBGKdynamics( T omega, Momenta<T,Lattice>& momenta, T singleScatAlbedo, T extinctionCoeff);
  /// Clone the object on its dynamic type.
  virtual AdDiAnisoBGKdynamics<T, Lattice>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium( int iPop, T rho2, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  virtual void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  virtual void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
                              LatticeStatistics<T>& statistics );
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega( T omega ) override;
private:
  T _omega;
  T _singleScatAlbedo;
  T _extinctionCoeff;
};



// ========= the BGK advection diffusion Stokes drag dynamics with a Smagorinsky turbulence model ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, template<typename U> class Lattice>
class SmagorinskyParticleAdvectionDiffusionBGKdynamics : public olb::AdvectionDiffusionBGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyParticleAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics );
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics );
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_, T dx_, T dt_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFacto_r, T rho, T pi[util::TensorVal<Lattice<T> >::n] );

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};

// ========= the BGK advection diffusion Stokes drag dynamics  ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, template<typename U> class Lattice>
class ParticleAdvectionDiffusionBGKdynamics : public olb::AdvectionDiffusionBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ParticleAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics );
private:
  T omega;  ///< relaxation parameter
};

} // namespace olb

#endif
