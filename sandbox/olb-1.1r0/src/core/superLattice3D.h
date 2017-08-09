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
 * The description of a 3D super lattice -- header file.
 */

#ifndef SUPER_LATTICE_3D_H
#define SUPER_LATTICE_3D_H

#include <vector>
#include "blockLattice3D.h"
#include "blockLatticeView3D.h"
#include "communication/communicator3D.h"
#include "postProcessing.h"
#include "serializer.h"
#include "communication/superStructure3D.h"


/// All OpenLB code is contained in this namespace.
namespace olb {




template<typename T, template<typename U> class Lattice> class SuperLatticeF3D;
template<typename T, typename S> class AnalyticalF3D;
template<typename T> class SuperGeometry3D;
template<typename T> class LoadBalancer;
template<typename T> class CuboidGeometry3D;
template <typename T> class IndicatorSphere3D;
template<typename T> class Communicator3D;
template<typename T, template<typename U> class Lattice> class BlockLatticeView3D;


/// A super lattice combines a number of block lattices that are ordered
/// in a cuboid geometry.
/** The communication between the block lattices is done by two
 * communicators. One (_commStream) is responsible to provide the data for
 * the streaming the other (_commBC) for the non-local boundary conditions.
 * To simplify the code structure ghost cells in an overlap of size
 * (_overlap) is introduced. It depends on the non-locality of the
 * boundary conditions but is at least one because of the streaming
 *
 * The algorithm is parallelized with mpi. The load balancer (_loadBalancer)
 * distributes the block lattices to processes.
 *
 * WARNING: For unstructured grids there is an interpolation needed
 * for the method buffer_outData in cuboidNeighbourhood which is not
 * yet implemented! Moreover this class needs to be changed so
 * that the number of times steps for the collision and streaming is
 * is dependent of the refinement level.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Lattice>
class SuperLattice3D : public SuperStructure3D<T>, public BufferSerializable {

private:
  /// Lattices with ghost cell layer of size overlap
  std::vector<BlockLattice3D<T, Lattice> > _extendedBlockLattices;
  /// View of the lattices without overlap
  std::vector<BlockLatticeView3D<T, Lattice> > _blockLattices;

  /// This communicator handles the communication for the streaming
  Communicator3D<T> _commStream;
  /// This communicator handles the communication for the postprocessors
  Communicator3D<T> _commBC;
  /// Specifies if there is communication for non local boundary conditions
  /// needed. It is automatically switched on if overlapBC >= 1 by
  /// calling the constructor. (default =false)
  bool _commBC_on;
  /// Statistic of the super structure
  LatticeStatistics<T> *_statistics;
  /// Specifies if there is statistic calculated. It is always
  /// needed for the ConstRhoBGK dynamics. (default =true)
  bool _statistics_on;

public:
  /// Construction of a super lattice
  SuperLattice3D(CuboidGeometry3D<T>& cGeometry,
                 LoadBalancer<T>& lb, int overlapBC = 0);

  SuperLattice3D(SuperGeometry3D<T>& superGeometry);
  ~SuperLattice3D();
  /// Read and write access to a block lattice
  BlockLattice3D<T, Lattice>& getExtendedBlockLattice(int locIC)
  {
    return _extendedBlockLattices[locIC];
  }
  ;
  /// Read only access to a block lattice
  BlockLattice3D<T, Lattice> const& getExtendedBlockLattice(int locIC) const
  {
    return _extendedBlockLattices[locIC];
  }
  ;
  /// Read and write access to a lattice (block lattice view, one
  /// without overlap).
  BlockLatticeView3D<T, Lattice>& getBlockLattice(int locIC)
  {
    return _blockLattices[locIC];
  }
  ;
  /// Read only access to a lattice
  BlockLatticeView3D<T, Lattice> const& getBlockLattice(int locIC) const
  {
    return _blockLattices[locIC];
  }
  ;

  /// Read and write access to the boundary communicator
  Communicator3D<T>& get_commBC()
  {
    return _commBC;
  }
  ;
  /// Read only access to the boundary communicator
  Communicator3D<T> const& get_commBC() const
  {
    return _commBC;
  }
  ;

  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics();
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const;

  /// Write access to lattice cells that returns false if
  /// (iX/iY/iZ) is not in any of the cuboids
  bool set(T iX, T iY, T iZ, Cell<T, Lattice> const& cell);
  /// Read only access to lattice cells that returns false if
  /// (iX/iY/iZ) is not in any of the cuboids
  bool get(T iX, T iY, T iZ, Cell<T, Lattice>& cell) const;
  /// Read only access to lattice cells over the cuboid number
  /// and local coordinates   WARNING!!! NO ERROR HANDLING IMPLEMENTED!!!
  Cell<T,Lattice> get(int iC, T locX, T locY, T locZ) const;
  Cell<T,Lattice> get(std::vector<int> latticeR) const;

  /// Write access to the memory of the data of the super structure
  virtual bool* operator() (int iCloc, int iX, int iY, int iZ, int iData)
  {
    return (bool*)&getExtendedBlockLattice(iCloc).get(iX+this->_overlap, iY+this->_overlap, iZ+this->_overlap)[iData];
  };
  /// Read only access to the dim of the data of the super structure
  virtual int getDataSize() const
  {
    return Lattice<T>::q;
  };
  /// Read only access to the data type dim of the data of the super structure
  virtual int getDataTypeSize() const
  {
    return sizeof(T);
  };

  /// Initialize all lattice cells to become ready for simulation
  void initialize();
  /// Defines the dynamics on a domain with a particular material number
  void defineDynamics(SuperGeometry3D<T>& sGeometry, int material,
                      Dynamics<T, Lattice>* dynamics);
  /// Defines rho on a domain with a particular material number
  void defineRho(SuperGeometry3D<T>& sGeometry, int material, AnalyticalF3D<T,T>& rho);
  /// Defines u on a domain with a particular material number
  void defineU(SuperGeometry3D<T>& sGeometry, int material, AnalyticalF3D<T,T>& u);
  /// Defines an external field on a domain with a particular material number
  void defineRhoU(SuperGeometry3D<T>& sGeometry, int material, AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u);
  // Defines a Population on a domain with a particular material number
  void definePopulations(SuperGeometry3D<T>& sGeometry, int material, AnalyticalF3D<T,T>& Pop);
  /// Defines an external field on a domain with a particular material number
  void defineExternalField(SuperGeometry3D<T>& sGeometry, int material,
                           int fieldBeginsAt, int sizeOfField, AnalyticalF3D<T,T>& field);
  /// Defines an external field on a IndicatorCircle domain
  void defineExternalField(SuperGeometry3D<T>& sGeometry, IndicatorSphere3D<T>& indicator,
                           int fieldBeginsAt, int sizeOfField, AnalyticalF3D<T,T>& field);
  /// Defines an external field on a domain with a particular material number
  void defineExternalField(SuperGeometry3D<T>& sGeometry, int material,
                           int fieldBeginsAt, int sizeOfField, SuperLatticeF3D<T,Lattice>& field);
  void setExternalParticleField(SuperGeometry3D<T>& sGeometry, AnalyticalF3D<T,T>& velocity,
                                ParticleIndicatorF3D<T,T>& sIndicator);
  /// Initializes by equilibrium on a domain with a particular material number
  void iniEquilibrium(SuperGeometry3D<T>& sGeometry, int material,
                      AnalyticalF3D<T,T>& rho , AnalyticalF3D<T,T>& u);

  /// Apply collision step to a rectangular domain
  void collide(T x0, T x1, T y0, T y1, T z0, T z1);
  /// Apply collision step to the whole domain
  void collide();
  /// TO BE DONE: Apply collision step to a rectangular domain,
  /// with fixed velocity
  // void staticCollide(T x0, T x1, T y0, T y1, T z0, T z1,
  //                  TensorField3D<T,3> const& u);
  /// TO BE DONE: Apply collision step to the whole domain,
  /// with fixed velocity
  // void staticCollide(TensorField3D<T,3> const& u);
  /// Apply streaming step to a rectangular domain
  void stream(T x0, T x1, T y0, T y1, T z0, T z1);
  /// Apply streaming step to the whole domain
  void stream();
  /// TO BE DONE: Apply first collision, then streaming step
  /// to a rectangular domain
  // void collideAndStream(T x0, T x1, T y0, T y1, T z0, T z1);
  /// Apply first collision, then streaming step
  /// to the whole domain
  void collideAndStream();
  /// Subtract a constant offset from the density within the whole domain
  void stripeOffDensityOffset (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T offset );
  /// Subtract a constant offset from the density within a rect. domain
  void stripeOffDensityOffset(T offset);
  /// Switches Statistics on (default on)
  void statisticsOn()
  {
    _statistics_on = true;
  };
  /// Switches Statistics off (default on). That speeds up
  /// the execution time.
  void statisticsOff()
  {
    _statistics_on = false;
  };

  /// Adds a coupling generator for one partner superLattice
  template<template<typename U> class Slattice>
  void addLatticeCoupling(SuperGeometry3D<T>& sGeometry, int material,
                          LatticeCouplingGenerator3D<T, Lattice> const& lcGen,
                          SuperLattice3D<T,Slattice>& partnerLattice );
  /// Executes coupling generator for one partner superLattice
  void executeCoupling();

  /// Number of data blocks for the serializable interface
  virtual std::size_t getNblock() const;
  /// Binary size for the serializer
  virtual std::size_t getSerializableSize() const;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  virtual bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode);

private:
  /// Resets and reduce the statistics
  void reset_statistics();
};

} // namespace olb

#endif
