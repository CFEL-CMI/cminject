/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2016 Robin Trunk
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic version of the collision, which modifies the particle
 *  distribution functions, by Orestis Malaspinas.
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

#include "core/postProcessing.h"
#include "core/blockLattice3D.h"

namespace olb {

/**
* This class interpolates missing f_i from values
* near the boundary to get a more stable outflow
* condition for the density. It is assumed that the
* next two cells in direction of this f_i have
* viable values.
*/
template<typename T, template<typename U> class Lattice>
class ConvectionBoundaryProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  ConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                int z1_, int discreteNormalX_,
                                int discreteNormalY_, int discreteNormalZ_);
  virtual int extent() const
  {
    return 0;
  }
  virtual int extent(int whichDirection) const
  {
    return 0;
  }
  virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
  virtual void processSubDomain ( BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_);
private:
  int interpolationPop[Lattice<T>::q];
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, template<typename U> class Lattice>
class ConvectionBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  ConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
                                         int z0_, int z1_, int discreteNormalX_,
                                         int discreteNormalY_, int discreteNormalZ_);
  virtual PostProcessor3D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};

/**
* This class copies missing values in the
* external field from the neighbour in normal direction.
* Therefore it is assumed this neighbour is a fluid cell.
*/
template<typename T, template<typename U> class Lattice>
class ExtFieldBoundaryProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  ExtFieldBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                              int z1_, int discreteNormalX_, int discreteNormalY_,
                              int discreteNormalZ_, int offset_);
  virtual int extent() const
  {
    return 0;
  }
  virtual int extent(int whichDirection) const
  {
    return 0;
  }
  virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
  virtual void processSubDomain ( BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_);
private:
  int x0, x1, y0, y1, z0, z1;
  int discreteNormalX, discreteNormalY, discreteNormalZ;
  int offset;
  bool par;
};

template<typename T, template<typename U> class Lattice>
class ExtFieldBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  ExtFieldBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                       int z1_, int discreteNormalX_, int discreteNormalY_,
                                       int discreteNormalZ_, int offset_);
  virtual PostProcessor3D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
  int discreteNormalX, discreteNormalY, discreteNormalZ;
  int offset;
};

/**
* This class resets some values of the distribution
* on the boundary that can have arbitrary values
* to be zero and thus ensures a correct computation
* of the density that is about to leave the domain.
*/
template<typename T, template<typename U> class Lattice>
class ZeroDistributionBoundaryProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  ZeroDistributionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                      int z1_, int discreteNormalX_, int discreteNormalY_,
                                      int discreteNormalZ_);
  virtual int extent() const
  {
    return 0;
  }
  virtual int extent(int whichDirection) const
  {
    return 0;
  }
  virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
  virtual void processSubDomain ( BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_);
private:
  int resetPop[Lattice<T>::q];
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, template<typename U> class Lattice>
class ZeroDistributionBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  ZeroDistributionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
      int z0_, int z1_, int discreteNormalX_,
      int discreteNormalY_, int discreteNormalZ_);
  virtual PostProcessor3D<T,Lattice>* generate() const;
  virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};
}

