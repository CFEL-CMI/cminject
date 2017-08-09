/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Orestis Malaspinas, Jonas Latt
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


#ifndef EXTENDED_FINITE_DIFFERENCE_BOUNDARY_2D_H
#define EXTENDED_FINITE_DIFFERENCE_BOUNDARY_2D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice2D.h"
#include "boundaryCondition2D.h"


namespace olb {

/**
* This class computes the finite difference approximation to LB boundary conditions
* on a flat wall in 2D with all the terms of the CE expansion. The details
* on the computations can be found in the thesis of
* Jonas Latt, "Hydrodynamic limit of lattice Boltzmann equations",
* University of Geneva, (2007).
*/
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class ExtendedStraightFdBoundaryPostProcessor2D : public LocalPostProcessor2D<T,Lattice> {
public:
  ExtendedStraightFdBoundaryPostProcessor2D(int x0_, int x1_, int y0_, int y1_);
  virtual int extent() const
  {
    return 1;
  }
  virtual int extent(int whichDirection) const
  {
    return 1;
  }
  virtual void process(BlockLattice2D<T,Lattice>& blockLattice);
  virtual void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ );
private:
  template<int deriveDirection>
  void interpolateGradients(BlockLattice2D<T,Lattice> const& blockLattice,
                            T velDeriv[Lattice<T>::d], int iX, int iY) const;
  template<int deriveDirection>
  void interpolateGradients (
    BlockLattice2D<T,Lattice> const& blockLattice,T& rhoDeriv,
    int iX, int iY ) const;
private:
  int x0, x1, y0, y1;
};

template<typename T, template<typename U> class Lattice, int direction, int orientation>
class ExtendedStraightFdBoundaryProcessorGenerator2D : public PostProcessorGenerator2D<T,Lattice> {
public:
  ExtendedStraightFdBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_);
  virtual PostProcessor2D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator2D<T,Lattice>*  clone() const;
};


////////// Factory function for Extended Finite Difference BC ///////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeBoundaryCondition2D<T,Lattice>*
createExtendedFdBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block);

template<typename T, template<typename U> class Lattice>
OnLatticeBoundaryCondition2D<T,Lattice>*
createExtendedFdBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block)
{
  return createExtendedFdBoundaryCondition2D<T,Lattice,BGKdynamics<T,Lattice> >(block);
}

}  // namespace olb

#endif
