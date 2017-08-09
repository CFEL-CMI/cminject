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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ADVECTION_DIFFUSION_MOMENTA_H
#define ADVECTION_DIFFUSION_MOMENTA_H

#include "advectionDiffusionLatticeDescriptors.h"

namespace olb {

/// Standard computation of velocity momenta in the bulk
template<typename T, template<typename U> class Lattice>
struct AdvectionDiffusionBulkMomenta : public Momenta<T,Lattice> {
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

namespace instances {

template<typename T, template<typename U> class Lattice>
AdvectionDiffusionBulkMomenta<T,Lattice>& getAdvectionDiffusionBulkMomenta();

}

}

#endif
