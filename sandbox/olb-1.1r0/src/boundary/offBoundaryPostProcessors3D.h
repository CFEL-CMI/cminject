/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_POST_PROCESSORS_3D_H
#define OFF_BOUNDARY_POST_PROCESSORS_3D_H

#include "core/postProcessing.h"
#include "core/blockLattice3D.h"

namespace olb {

/**
* This class computes the Linear Bouzidi BC
*/

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBouzidiLinearPostProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  ZeroVelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
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
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
private:
  int x, y, z;
  int xN, yN, zN, xB, yB, zB;
  int iPop, opp, iPop2;
  T q, dist;
};

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBounceBackPostProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  ZeroVelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
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
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
private:
  int x, y, z;
  int xN, yN, zN;
  int iPop, opp;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class VelocityBouzidiLinearPostProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  VelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
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
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
private:
  int x, y, z;
  int xN, yN, zN, xB, yB, zB;
  int iPop, opp, iPop2;
  T q, dist;
  T ufrac;
};

template<typename T, template<typename U> class Lattice>
class VelocityBounceBackPostProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  VelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
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
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
private:
  int x, y, z;
  int xN, yN, zN;
  int iPop, opp;
  T dist;
};

/**
* Linear Bouzidi BC Generator
*/

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBouzidiLinearPostProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  ZeroVelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  virtual PostProcessor3D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
  int x, y, z;
  int iPop;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBounceBackPostProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  ZeroVelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  virtual PostProcessor3D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
  int x, y, z;
  int iPop;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class VelocityBouzidiLinearPostProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  VelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  virtual PostProcessor3D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
  int x, y, z;
  int iPop;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class VelocityBounceBackPostProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  VelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  virtual PostProcessor3D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
  int x, y, z;
  int iPop;
  T dist;
};

}

#endif
