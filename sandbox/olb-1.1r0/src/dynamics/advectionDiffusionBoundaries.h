/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_H
#define ADVECTION_DIFFUSION_BOUNDARIES_H


#include "advectionDiffusionLatticeDescriptors.h"
#include "advectionDiffusionDynamics.h"
#include "dynamics/dynamics.h"

namespace olb {

/**
* This class computes the Advection Diffusion BC with general dynamics.
*/
//===================================================================================
//================= AdvectionDiffusionDynamcison Flat Boundaries =========
//===================================================================================

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
class AdvectionDiffusionBoundariesDynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  AdvectionDiffusionBoundariesDynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual AdvectionDiffusionBoundariesDynamics<T, Lattice, Dynamics, direction, orientation>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
private:
  Dynamics boundaryDynamics;
};

//===================================================================================
//================= AdvectionDiffusionDynamcis On Edges =========
//===================================================================================

template<typename T, template<typename U> class Lattice, typename Dynamics, int plane, int normal1, int normal2>
class AdvectionDiffusionEdgesDynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  AdvectionDiffusionEdgesDynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual AdvectionDiffusionEdgesDynamics<T, Lattice, Dynamics, plane, normal1, normal2>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
private:
  Dynamics boundaryDynamics;
};


//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 2D Boundaries =========
//===================================================================================

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal>
class AdvectionDiffusionCornerDynamics2D : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  AdvectionDiffusionCornerDynamics2D(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual AdvectionDiffusionCornerDynamics2D<T, Lattice, Dynamics, xNormal, yNormal>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
private:
  Dynamics boundaryDynamics;
};

//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 3D Boundaries =========
//===================================================================================

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal, int zNormal>
class AdvectionDiffusionCornerDynamics3D : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  AdvectionDiffusionCornerDynamics3D(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual AdvectionDiffusionCornerDynamics3D<T, Lattice, Dynamics, xNormal, yNormal, zNormal>* clone() const;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
private:
  Dynamics boundaryDynamics;
};



}  // namespace olb

#endif
