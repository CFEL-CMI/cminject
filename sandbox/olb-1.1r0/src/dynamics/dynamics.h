/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Mathias J. Krause
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
#ifndef LB_DYNAMICS_H
#define LB_DYNAMICS_H

#include "latticeDescriptors.h"
#include "core/util.h"
#include "core/postProcessing.h"
#include "core/latticeStatistics.h"

namespace olb {

namespace dynamicParams {
// Use 0-99 for relaxation parameters
const int omega_shear = 0;
const int omega_bulk  = 1;

// Use 100-199 for material constants
const int sqrSpeedOfSound = 100; // Speed of sound squared
const int sqrInvSpeedOfSound = 101; // Inverse speed of sound squared

// Use 1000 and higher for custom user-defined constants
}

template<typename T, template<typename U> class Lattice> class Cell;

/// Interface for the dynamics classes
template<typename T, template<typename U> class Lattice>
struct Dynamics {
  /// Destructor: virtual to enable inheritance
  virtual ~Dynamics() { }
  /// Clone the object on its dynamic type.
  virtual Dynamics<T,Lattice>* clone() const =0;
  /// Implementation of the collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) =0;
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) =0;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const =0;
  /// Initialize cell at equilibrium distribution
  virtual void iniEquilibrium(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d]);
  /// Compute particle density on the cell.
  /** \return particle density
   */
  virtual T computeRho(Cell<T,Lattice> const& cell) const =0;
  /// Compute fluid velocity on the cell.
  /** \param u fluid velocity
   */
  virtual void computeU( Cell<T,Lattice> const& cell,
                         T u[Lattice<T>::d] ) const =0;
  /// Compute fluid momentum (j=rho*u) on the cell.
  /** \param j fluid momentum
   */
  virtual void computeJ( Cell<T,Lattice> const& cell,
                         T j[Lattice<T>::d] ) const =0;
  /// Compute components of the stress tensor on the cell.
  /** \param pi stress tensor */
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const =0;
  /// Compute fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const =0;
  /// Compute all momenta on the cell, up to second order.
  /** \param rho particle density
   *  \param u fluid velocity
   *  \param pi stress tensor
   */
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const =0;
  /// Access particle populations through the dynamics object.
  /** Default implementation: access cell directly.
   */
  virtual void computePopulations(Cell<T,Lattice> const& cell, T* f) const;
  /// Access external fields through the dynamics object.
  /** Default implementation: access cell directly.
   */
  virtual void computeExternalField (
    Cell<T,Lattice> const& cell, int pos, int size, T* ext ) const;
  /// Set particle density on the cell.
  /** \param rho particle density
   */
  virtual void defineRho(Cell<T,Lattice>& cell, T rho) =0;
  virtual void defineRho(int iPop, T rho);
  /// Set fluid velocity on the cell.
  /** \param u fluid velocity
   */
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) =0;
  /// Functions for offLattice Velocity boundary conditions
  virtual void setBoundaryIntersection(int iPop, T distance);
  virtual bool getBoundaryIntersection(int iPop, T point[Lattice<T>::d]);
  virtual void defineU(const T u[Lattice<T>::d]);
  virtual void defineU(int iPop, const T u[Lattice<T>::d]);
  virtual T    getVelocityCoefficient(int iPop);

  /// Define fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) =0;
  /// Define all momenta on the cell, up to second order.
  /** \param rho particle density
   *  \param u fluid velocity
   *  \param pi stress tensor
   */
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) =0;
  /// Define particle populations through the dynamics object.
  /** Default implementation: access cell directly.
   */
  virtual void definePopulations(Cell<T,Lattice>& cell, const T* f);
  /// Define external fields through the dynamics object.
  /** Default implementation: access cell directly.
   */
  virtual void defineExternalField(Cell<T,Lattice>& cell, int pos, int size, const T* ext);
  /// Add external fields through the dynamics object.
  /** Similar to defineExternalField(),but instead of replacing existing values
   *  ext is added to existing values.
   */
  virtual void addExternalField(Cell<T,Lattice>& cell, int pos, int size, const T* ext);
  /// Add external fields through the dynamics object.
  /** Similar to defineExternalField(),but instead of replacing existing values
   *  ext is multiplied to existing values.
   */
  virtual void multiplyExternalField(Cell<T,Lattice>& cell, int pos, int size, const T* ext);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const =0;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega) =0;
  /// Get local value of any parameter
  virtual T getParameter(int whichParameter) const;
  /// Set local value of any parameter
  virtual void setParameter(int whichParameter, T value);
};

/// Interface for classes that compute velocity momenta
/** This class is useful for example to distinguish between bulk and
 * boundary nodes, given that on the boundaries, a particular strategy
 * must be applied to compute the velocity momenta.
 */
template<typename T, template<typename U> class Lattice>
struct Momenta {
  /// Destructor: virtual to enable inheritance
  virtual ~Momenta() { }
  /// Compute particle density on the cell.
  virtual T computeRho(Cell<T,Lattice> const& cell) const =0;
  /// Compute fluid velocity on the cell.
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const =0;
  /// Compute fluid momentum on the cell.
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const =0;
  /// Compute components of the stress tensor on the cell.
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const =0;
  /// Compute fluid velocity and particle density on the cell.
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  /// Compute all momenta on the cell, up to second order.
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Set particle density on the cell.
  virtual void defineRho(Cell<T,Lattice>& cell, T rho) =0;
  /// Set fluid velocity on the cell.
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) =0;
  /// Define fluid velocity and particle density on the cell.
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Define all momenta on the cell, up to second order.
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) =0;
};

/// Abstract base for dynamics classes
/** In this version of the Dynamics classes, the computation of the
 * velocity momenta is taken care of by an object of type Momenta.
 */
template<typename T, template<typename U> class Lattice>
class BasicDynamics : public Dynamics<T,Lattice> {
public:
  /// Must be contructed with an object of type Momenta
  BasicDynamics(Momenta<T,Lattice>& momenta);
  /// Implemented via the Momenta object
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Implemented via the Momenta object
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Implemented via the Momenta object
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const;
  /// Implemented via the Momenta object
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Implemented via the Momenta object
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  /// Implemented via the Momenta object
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Implemented via the Momenta object
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Implemented via the Momenta object
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]);
  /// Implemented via the Momenta object
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Implemented via the Momenta object
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] );
protected:
  Momenta<T,Lattice>& _momenta;  ///< computation of velocity momenta
};

/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class BGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  BGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Clone the object on its dynamic type.
  virtual BGKdynamics<T,Lattice>* clone() const;
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
  virtual void setOmega(T omega);
private:
  T _omega;  ///< relaxation parameter
};

/// Implementation of the pressure-corrected BGK collision step
template<typename T, template<typename U> class Lattice>
class ConstRhoBGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  ConstRhoBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Clone the object on its dynamic type.
  virtual ConstRhoBGKdynamics<T,Lattice>* clone() const;
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
  virtual void setOmega(T omega);
private:
  T _omega;  ///< relaxation parameter
};

/// Implementation of the so-called incompressible collision step
template<typename T, template<typename U> class Lattice>
class IncBGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  IncBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Clone the object on its dynamic type.
  virtual IncBGKdynamics<T,Lattice>* clone() const;
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
  virtual void setOmega(T omega);
private:
  T _omega;  ///< relaxation parameter
};



/// Implementation of the Regularized BGK collision step
/** This model is substantially more stable than plain BGK, and has roughly
 * the same efficiency. However, it cuts out the modes at higher Knudsen
 * numbers and can not be used in the regime of rarefied gases.
 */
template<typename T, template<typename U> class Lattice>
class RLBdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  RLBdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Clone the object on its dynamic type.
  virtual RLBdynamics<T,Lattice>* clone() const;
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
  virtual void setOmega(T omega);
private:
  T _omega;  ///< relaxation parameter
};

/// Implementation of Regularized BGK collision, followed by any Dynamics
template<typename T, template<typename U> class Lattice, typename Dynamics>
class CombinedRLBdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  CombinedRLBdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Clone the object on its dynamic type.
  virtual CombinedRLBdynamics<T, Lattice, Dynamics>* clone() const;
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
  virtual void setOmega(T omega);
  /// Get local value of any parameter
  virtual T getParameter(int whichParameter) const;
  /// Set local value of any parameter
  virtual void setParameter(int whichParameter, T value);
private:
  Dynamics _boundaryDynamics;
};

/// Implementation of the BGK collision step with external force
template<typename T, template<typename U> class Lattice>
class ForcedBGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  ForcedBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Clone the object on its dynamic type.
  virtual ForcedBGKdynamics<T,Lattice>* clone() const;
  ///  Compute fluid velocity on the cell.
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Compute fluid velocity and particle density on the cell.
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
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
  virtual void setOmega(T omega);
protected:
  T _omega;  ///< relaxation parameter
  static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
  static const int sizeOfForce   = Lattice<T>::ExternalField::sizeOfForce;
};

/// Implementation of the BGK collision step with external force
template<typename T, template<typename U> class Lattice>
class ResettingForcedBGKdynamics : public ForcedBGKdynamics<T,Lattice> {
public:
  ResettingForcedBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  inline void setForce(T force[3])
  {
//    _frc[0] = force[0];
//    _frc[1] = force[1];
//    _frc[2] = force[2];
    _frc[0] = 0.0;
    _frc[1] = 0.0;
    _frc[2] = 0.0;
  }
private:
  T _frc[3];
};

/// Other Implementation of the BGK collision step with external force
template<typename T, template<typename U> class Lattice>
class ForcedShanChenBGKdynamics : public ForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ForcedShanChenBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  ///  Compute fluid velocity on the cell.
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Compute fluid velocity and particle density on the cell.
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
};

/// Implementation of the 3D D3Q13 dynamics
/** This is (so far) the minimal existing 3D model, with only 13
 * directions. Three different relaxation times are used to achieve
 * asymptotic hydrodynamics, isotropy and galilean invariance.
 */
template<typename T, template<typename U> class Lattice>
class D3Q13dynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  D3Q13dynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Clone the object on its dynamic type.
  virtual D3Q13dynamics<T,Lattice>* clone() const;
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
  virtual void setOmega(T omega);
private:
  T lambda_nu;        ///< first relaxation parameter
  T lambda_nu_prime;  ///< second relaxation parameter
};

/// Standard computation of velocity momenta in the bulk
template<typename T, template<typename U> class Lattice>
struct BulkMomenta : public Momenta<T,Lattice> {
  /// Compute particle density on the cell.
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Compute fluid velocity on the cell.
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Compute fluid momentum on the cell.
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const;
  /// Compute components of the stress tensor on the cell.
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Compute fluid velocity and particle density on the cell.
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  /// Compute all momenta on the cell, up to second order.
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Set particle density on the cell.
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Set fluid velocity on the cell.
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]);
  /// Define fluid velocity and particle density on the cell.
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Define all momenta on the cell, up to second order.
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] );
};

/// Velocity is stored in external scalar (and computed e.g. in a PostProcessor)
template<typename T, template<typename U> class Lattice>
struct ExternalVelocityMomenta : public Momenta<T,Lattice> {
  /// Compute particle density on the cell.
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Compute fluid velocity on the cell.
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Compute fluid momentum on the cell.
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const;
  /// Compute components of the stress tensor on the cell.
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Compute fluid velocity and particle density on the cell.
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  /// Compute all momenta on the cell, up to second order.
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Set particle density on the cell.
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Set fluid velocity on the cell.
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]);
  /// Define fluid velocity and particle density on the cell.
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Define all momenta on the cell, up to second order.
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] );
};

/// Implementation of "bounce-back" dynamics
/** This is a very popular way to implement no-slip boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Lattice>
class BounceBack : public Dynamics<T,Lattice> {
public:
  /// A fictitious density value on bounce-back in not fixed on nodes via this constructor.
  BounceBack();
  /// You may fix a fictitious density value on bounce-back nodes via this constructor.
  BounceBack(T rho);
  /// Clone the object on its dynamic type.
  virtual BounceBack<T,Lattice>* clone() const;
  /// Yields 0;
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Yields 1;
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Yields 0;
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Yields 0;
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const;
  /// Yields NaN
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Does nothing
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Does nothing
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]);
  /// Does nothing
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Does nothing
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] );
  /// Yields NaN
  virtual T getOmega() const;
  /// Does nothing
  virtual void setOmega(T omega);
private:
  T _rho;
  bool _rhoFixed;
};

/// Implementation of "bounce-back velocity" dynamics
/** This is a very popular way to implement no-slip boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics. It
 * fixes the velociy to a given velocity _u.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Lattice>
class BounceBackVelocity : public Dynamics<T,Lattice> {
public:
  /// A fictitious density value on bounce-back in not fixed on nodes via this constructor.
  BounceBackVelocity(const T u[Lattice<T>::d]);
  /// You may fix a fictitious density value on bounce-back nodes via this constructor.
  BounceBackVelocity(const T rho, const T u[Lattice<T>::d]);
  /// Clone the object on its dynamic type.
  virtual BounceBackVelocity<T,Lattice>* clone() const;
  /// Yields 0;
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step, bounce back with a fixed velocity _u
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Retuns rho (if defined else zero)
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Retuns _u
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Retuns rho (if defined else zero) times _u
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const;
  /// Yields NaN
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Retuns rho (if defined else zero) and _u
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Devines the velocity rho
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Devines the velocity _u
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]);
  /// Devines rho and _u
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Does nothing
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] );
  /// Yields NaN
  virtual T getOmega() const;
  /// Does nothing
  virtual void setOmega(T omega);
private:
  T _rho;
  bool _rhoFixed;
  T _u[Lattice<T>::d];
};

/// Implementation of "bounce-back anti" dynamics
/** This is a way to implement a Dirichlet rho/pressure boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics. It
 * fixes the rho to a given _rho.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Lattice>
class BounceBackAnti : public Dynamics<T,Lattice> {
public:
  /// A fictitious density value on bounce-back in not fixed on nodes via this constructor.
  BounceBackAnti();
  /// You may fix a fictitious density value on bounce-back nodes via this constructor.
  BounceBackAnti(T rho);
  /// Clone the object on its dynamic type.
  virtual BounceBackAnti<T,Lattice>* clone() const;
  /// Yields 0;
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step, bounce back with a fixed velocity _u
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Retuns rho (if defined else zero)
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Retuns _u
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Retuns rho (if defined else zero) times _u
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const;
  /// Yields NaN
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Retuns rho (if defined else zero) and _u
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Devines the velocity rho
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Devines the velocity _u
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]);
  /// Devines rho and _u
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Does nothing
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] );
  /// Yields NaN
  virtual T getOmega() const;
  /// Does nothing
  virtual void setOmega(T omega);
private:
  T _rho;
  bool _rhoFixed;
  T _u[Lattice<T>::d];
};

/** Robin Boundary for Diffusion Equation
 * documentation Hiorth, Lad, Evje and Skjæveland 2008
 */
template<typename T, template<typename U> class Lattice>
class BounceBackReflective : public BounceBack<T,Lattice> {
public:
  BounceBackReflective(T zeta);
  virtual T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) final;
private:
  T _zeta;
};


/// Implementation of a "dead cell" that does nothing
template<typename T, template<typename U> class Lattice>
class NoDynamics : public Dynamics<T,Lattice> {
public:
  /// You may fix a fictitious density value on no dynamics node via this constructor.
  NoDynamics(T rho = T(1) );
  /// Clone the object on its dynamic type.
  virtual NoDynamics<T,Lattice>* clone() const;
  /// Yields 0;
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Yields 1;
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Yields 0;
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const;
  /// Yields 0;
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const;
  /// Yields NaN
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Does nothing
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Does nothing
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]);
  /// Does nothing
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Does nothing
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] );
  /// Yields NaN
  virtual T getOmega() const;
  /// Does nothing
  virtual void setOmega(T omega);

private:
  /// Default rho=1
  T _rho;
};

/// Dynamics for offLattice boundary conditions
/// OffDynamics are basically NoDynamics with the additional functionality
/// to store given velocities exactly at boundary links.
template<typename T, template<typename U> class Lattice>
class OffDynamics : public NoDynamics<T,Lattice> {
public:
  /// Constructor
  OffDynamics(const T _location[Lattice<T>::d]);
  /// Constructor
  OffDynamics(const T _location[Lattice<T>::d], T _distances[Lattice<T>::q]);
  /// Returns local stored rho which is updated if the bc is used as velocity!=0 condition
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  /// Returns an average of the locally stored u
  virtual void computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const;
  /// Set Intersection of the link and the boundary
  virtual void setBoundaryIntersection(int iPop, T distance);
  /// Get Intersection of the link and the boundary
  virtual bool getBoundaryIntersection(int iPop, T intersection[Lattice<T>::d]);
  /// Set particle density on the cell.
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  /// Set single velocity
  virtual void defineRho(int iPop, T rho);
  /// Set fluid velocity on the cell.
  virtual void defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]);
  /// Set constant velocity
  virtual void defineU(const T u[Lattice<T>::d]);
  /// Set single velocity
  virtual void defineU(int iPop, const T u[Lattice<T>::d]);
  /// Get VelocitySummand for Bouzidi-Boundary Condition
  virtual T getVelocityCoefficient(int iPop);

private:
  T _rho;
  T _u[Lattice<T>::q][Lattice<T>::d];
  T location[Lattice<T>::d];
  T distances[Lattice<T>::q];
  T boundaryIntersection[Lattice<T>::q][Lattice<T>::d];
  T velocityCoefficient[Lattice<T>::q];
};

/// Implementation of density sink by setting a zero distribution on the cell
template<typename T, template<typename U> class Lattice>
class ZeroDistributionDynamics : public NoDynamics<T,Lattice> {
public:
  /// Constructor.
  ZeroDistributionDynamics();
  /// Clone the object on its dynamic type.
  virtual ZeroDistributionDynamics<T,Lattice>* clone() const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Yields 1
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
};


namespace instances {

template<typename T, template<typename U> class Lattice>
BulkMomenta<T,Lattice>& getBulkMomenta();

template<typename T, template<typename U> class Lattice>
ExternalVelocityMomenta<T,Lattice>& getExternalVelocityMomenta();

template<typename T, template<typename U> class Lattice>
BounceBack<T,Lattice>& getBounceBack();

template<typename T, template<typename U> class Lattice>
BounceBackVelocity<T,Lattice>& getBounceBackVelocity(const double rho, const double u[Lattice<T>::d]);

template<typename T, template<typename U> class Lattice>
BounceBackAnti<T,Lattice>& getBounceBackAnti(const double rho);

template<typename T, template<typename U> class Lattice>
NoDynamics<T,Lattice>& getNoDynamics(T rho = T(1) );

template<typename T, template<typename U> class Lattice>
ZeroDistributionDynamics<T,Lattice>& getZeroDistributionDynamics();

}

}

#endif
