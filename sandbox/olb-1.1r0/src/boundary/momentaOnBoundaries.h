/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Implementation of boundary cell dynamics -- header file.
 */
#ifndef MOMENTA_ON_BOUNDARIES_H
#define MOMENTA_ON_BOUNDARIES_H

#include "dynamics/dynamics.h"
#include "core/util.h"
#include "core/cell.h"

namespace olb {

template<typename T, template<typename U> class Lattice> class Cell;

/// Dirichlet condition on velocity and/or pressure
template<typename T, template<typename U> class Lattice>
class DirichletBoundaryMomenta : public Momenta<T,Lattice> {
};

template<typename T, template<typename U> class Lattice>
class EquilibriumBM : public DirichletBoundaryMomenta<T,Lattice> {
public:
  EquilibriumBM();
  EquilibriumBM( T rho, const T u[Lattice<T>::d] );
  virtual T computeRho( Cell<T,Lattice> const& cell ) const;
  virtual void computeU( Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const;
  virtual void computeJ( Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const;
  virtual void computeStress( Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
                              T pi[util::TensorVal<Lattice<T> >::n] ) const;
  virtual void defineRho( Cell<T,Lattice>& cell, T rho);
  virtual void defineU( Cell<T,Lattice>& cell, const T u[Lattice<T>::d]);
  virtual void defineAllMomenta( Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
                                 const T pi[util::TensorVal<Lattice<T> >::n] );
private:
  T _rho;
  T _u[Lattice<T>::d];   ///< value of the velocity on the boundary
};

/// Computation of velocity momenta on a velocity boundary
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class VelocityBM : virtual public DirichletBoundaryMomenta<T,Lattice> {
public:
  /// Default Constructor: initialization to zero
  VelocityBM();
  /// Constructor with boundary initialization
  VelocityBM(const T u[Lattice<T>::d]);

  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  virtual void computeU ( Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const;
  virtual void computeJ ( Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const;
  void computeU(T u[Lattice<T>::d]) const;
  virtual void defineRho(Cell<T,Lattice>& cell, T rho) ;
  virtual void defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]) ;
  void defineU(const T u[Lattice<T>::d]);
  virtual void defineAllMomenta( Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
                                 const T pi[util::TensorVal<Lattice<T> >::n] );
private:
  T _u[Lattice<T>::d];   ///< value of the velocity on the boundary
};

/// Computation of velocity momenta on a velocity boundary
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class PressureBM : virtual public DirichletBoundaryMomenta<T,Lattice> {
public:
  /// Default Constructor: initialization to u=0, rho=1
  PressureBM();
  /// Constructor with boundary initialization
  PressureBM(const T values_[Lattice<T>::d]);

  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  T computeRho() const;
  virtual void computeU( Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const;
  virtual void computeJ( Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const;
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  void defineRho(T rho);
  virtual void defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]);
  virtual void defineAllMomenta( Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
                                 const T pi[util::TensorVal<Lattice<T> >::n] );
private:
  /// Velocity/Density on boundary.
  /** Contains velocity on the boundary, except for values[direction] that
   * contains a prescription for the density.
   */
  T _values[Lattice<T>::d];
};

/// Here, the stress is computed from the particle distribution functions
template<typename T, template<typename U> class Lattice>
class FreeStressBM : virtual public DirichletBoundaryMomenta<T,Lattice> {
public:
  virtual void computeStress( Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
                              T pi[util::TensorVal<Lattice<T> >::n] ) const;
};

/// Use special trick to compute u resp. rho, but compute pi from part. distr. functions
template<typename T, template<typename U> class Lattice,
         template <
           typename T_, template<typename U_> class Lattice_,
           int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
class BasicDirichletBM : public FreeStressBM<T, Lattice>, public HydroBM<T,Lattice,direction,orientation> {
};

/// Computation of the stress tensor for regularized boundary
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class RegularizedBM : virtual public DirichletBoundaryMomenta<T,Lattice> {
public:
  /// Stress tensor
  virtual void computeStress( Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
                              T pi[util::TensorVal<Lattice<T> >::n] ) const;
};

/// Regularized velocity boundary node
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class RegularizedVelocityBM : public RegularizedBM<T,Lattice,direction,orientation>, public VelocityBM<T,Lattice,direction,orientation> {
public:
  RegularizedVelocityBM() { }
  RegularizedVelocityBM(const T u[Lattice<T>::d]) : VelocityBM<T,Lattice,direction,orientation>(u) {}
};

/// Regularized pressure boundary node
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class RegularizedPressureBM : public RegularizedBM<T,Lattice,direction,orientation>, public PressureBM<T,Lattice,direction,orientation> {
public:
  RegularizedPressureBM() { }
  RegularizedPressureBM(const T values[Lattice<T>::d]) : PressureBM<T,Lattice,direction,orientation>(values) { }
};

/// In this class, the velocity is fixed
/**
 * As opposed to VelocityBM, the pressure is however not
 * computed from a special trick on the boundary, but the
 * same way it would be in the bulk.
 */
template<typename T, template<typename U> class Lattice>
class FixedVelocityBM : public Momenta<T,Lattice> {
public:
  virtual T computeRho(Cell<T,Lattice> const& cell) const;
  virtual void computeU( Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const;
  virtual void computeJ( Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const;
  virtual void computeStress( Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
                              T pi[util::TensorVal<Lattice<T> >::n] ) const;
  virtual void computeRhoU( Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d]) const;
  virtual void computeAllMomenta( Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d],
                                  T pi[util::TensorVal<Lattice<T> >::n] ) const;
  virtual void defineRho(Cell<T,Lattice>& cell, T rho);
  virtual void defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]);
  virtual void defineRhoU( Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d]);
  virtual void defineAllMomenta( Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
                                 const T pi[util::TensorVal<Lattice<T> >::n] );
private:
  BulkMomenta<T,Lattice> _basicMomenta;
  T _fixU[Lattice<T>::d];
};

}  // namespace olb


#endif

