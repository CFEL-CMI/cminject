/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
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

#ifndef NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_H
#define NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_H

#include <cmath>

#include "core/spatiallyExtendedObject3D.h"
#include "core/postProcessing.h"
#include "core/blockLattice3D.h"
#include "core/units.h"
#include "advectionDiffusionForces.hh"
#include "advectionDiffusionForces.h"

namespace olb {

/**
* Multiphysics class for coupling between different lattices.
*/
//======================================================================
// ========Regularized NSDiffusion Coupling 3D ====================//
//======================================================================
template<typename T, template<typename U> class Lattice>
class NavierStokesAdvectionDiffusionCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,Lattice> {
public:
  NavierStokesAdvectionDiffusionCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<SpatiallyExtendedObject3D* > partners_);
  virtual int extent() const
  {
    return 1;
  }
  virtual int extent(int whichDirection) const
  {
    return 1;
  }
  virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
  virtual void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_);
private:
  int x0, x1, y0, y1, z0, z1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;

  std::vector<SpatiallyExtendedObject3D*> partners;
};

template<typename T, template<typename U> class Lattice>
class NavierStokesAdvectionDiffusionCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,Lattice> {
public:
  NavierStokesAdvectionDiffusionCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  virtual PostProcessor3D<T,Lattice>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const;
  virtual LatticeCouplingGenerator3D<T,Lattice>* clone() const;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};

//==================================================================================================
// ========Coupling 3D of Navier-Stokes on Advection-Diffusion with Stokes drag====================//
//==================================================================================================
template<typename T, template<typename U> class Lattice>
class AdvectionDiffusionParticleCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,Lattice> {
public:
  AdvectionDiffusionParticleCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_, int offset_,
    std::vector<SpatiallyExtendedObject3D* > partners_, std::vector<std::reference_wrapper<advectionDiffusionForce3D<T, Lattice> > > forces_);
  virtual int extent() const
  {
    return 1;
  }
  virtual int extent(int whichDirection) const
  {
    return 1;
  }
  virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
  virtual void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_);
private:
  int x0, x1, y0, y1, z0, z1, iC;
  T dragCoeff;
  int offset;
  T *vel, *vel_new, *velXp, *velXn, *velYp, *velYn, *velZp, *velZn;
  bool par = true;
  Cell<T,descriptors::particleAdvectionDiffusionD3Q7Descriptor> *adCell;
  Cell<T,Lattice> *nsCell;

protected:
  std::vector<std::reference_wrapper<advectionDiffusionForce3D<T, Lattice> > > forces;
};

template<typename T, template<typename U> class Lattice>
class AdvectionDiffusionParticleCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,Lattice> {
public:
  AdvectionDiffusionParticleCouplingGenerator3D(int offset_);
  virtual PostProcessor3D<T,Lattice>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const;
  virtual LatticeCouplingGenerator3D<T,Lattice>* clone() const;
  void addForce(advectionDiffusionForce3D<T,Lattice> &force);

private:
  int offset;

protected:
  std::vector<std::reference_wrapper<advectionDiffusionForce3D<T, Lattice> > > ADforces;
};

}

#endif
