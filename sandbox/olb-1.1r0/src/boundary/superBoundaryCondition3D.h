/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
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
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef SUPER_BOUNDARY_CONDITION_3D_H
#define SUPER_BOUNDARY_CONDITION_3D_H

#include <vector>
#include "boundaryCondition3D.h"
#include "dynamics/advectionDiffusionBoundaryCondition3D.h"
#include "geometry/superGeometryStatistics3D.h"
#include "core/superLattice3D.h"
#include "io/ostreamManager.h"
#include "extendedFiniteDifferenceBoundary3D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

/// A helper for initialising 3D boundaries for super lattices.
/** Here we have methods that initializes the local postprocessors and the
 * communicator (_commBC in SuperLattice) for boundary conditions
 * for a given global point or global range.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Lattice>
class sOnLatticeBoundaryCondition3D {

public:
  /// Constructor
  sOnLatticeBoundaryCondition3D(SuperLattice3D<T, Lattice>& sLattice);
  /// Copy construction
  sOnLatticeBoundaryCondition3D(
    sOnLatticeBoundaryCondition3D<T, Lattice> const& rhs);
  /// Copy assignment
  sOnLatticeBoundaryCondition3D operator=(sOnLatticeBoundaryCondition3D<T,
                                          Lattice> rhs);
  /// Destructor
  ~sOnLatticeBoundaryCondition3D();

  void addVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);
  void addSlipBoundary(SuperGeometry3D<T>& superGeometry, int material);
  void addPressureBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);
  void addConvectionBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega, T* uAv=NULL);
  void addTemperatureBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);
  void addDiffuseReflectionBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega, T zeta);
  void addConvectionBoundary(SuperGeometry3D<T>& superGeometry, int material);
  void addExtFieldBoundary(SuperGeometry3D<T>& superGeometry, int material, int offset);
  void addZeroDistributionBoundary(SuperGeometry3D<T>& superGeometry, int material);

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material);

  SuperLattice3D<T, Lattice>& getSuperLattice()
  {
    return _sLattice;
  };

  std::vector<OnLatticeBoundaryCondition3D<T, Lattice>*>& getBlockBCs()
  {
    return _blockBCs;
  };

  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Lattice>*>& getADblockBCs()
  {
    return _ADblockBCs;
  };

  int getOverlap()
  {
    return _overlap;
  };

  void setOverlap(int overlap)
  {
    _overlap = overlap;
  };

  void outputOn();
  void outputOff();


private:
  mutable OstreamManager clout;
  SuperLattice3D<T, Lattice>& _sLattice;
  std::vector<OnLatticeBoundaryCondition3D<T, Lattice>*> _blockBCs;
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Lattice>*> _ADblockBCs;
  int _overlap;
  bool _output;
};

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createLocalBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC);

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createInterpBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC);

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createExtFdBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC);


template<typename T, template<typename U> class Lattice>
void createLocalBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC)
{
  createLocalBoundaryCondition3D<T, Lattice, RLBdynamics<T, Lattice> > (sBC);
}
template<typename T, template<typename U> class Lattice>
void createInterpBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC)
{
  createInterpBoundaryCondition3D<T, Lattice, BGKdynamics<T, Lattice> > (sBC);
}
template<typename T, template<typename U> class Lattice>
void createExtFdBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC)
{
  createInterpBoundaryCondition3D<T, Lattice, BGKdynamics<T, Lattice> > (sBC);
}

} // namespace olb

#endif
