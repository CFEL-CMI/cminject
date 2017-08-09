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

/** \file A helper for initialising 3D boundaries -- header file.  */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H

#include "advectionDiffusionBoundaryCondition3D.h"
//#include "advectionDiffusionBoundaryCondition3D.hh"
//#include "advectionDiffusionBoundaryPostProcessor3D.hh"

namespace olb {

template<typename T, template<typename U> class Lattice, class BoundaryManager>
class AdvectionDiffusionBoundaryConditionInstantiator3D : public OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Lattice> {
public:
  AdvectionDiffusionBoundaryConditionInstantiator3D( BlockLatticeStructure3D<T,Lattice>& block_ );
  ~AdvectionDiffusionBoundaryConditionInstantiator3D();

  void addTemperatureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega);

  // Temperature Boundary Conditions for edges ...
  void addTemperatureBoundaryEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  void addTemperatureBoundaryEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega);

  // Temperature Boundary Conditions for Corners ...
  void addTemperatureBoundaryCornerNNN(int x, int y, int z, T omega);
  void addTemperatureBoundaryCornerNNP(int x, int y, int z, T omega);
  void addTemperatureBoundaryCornerNPN(int x, int y, int z, T omega);
  void addTemperatureBoundaryCornerNPP(int x, int y, int z, T omega);
  void addTemperatureBoundaryCornerPNN(int x, int y, int z, T omega);
  void addTemperatureBoundaryCornerPNP(int x, int y, int z, T omega);
  void addTemperatureBoundaryCornerPPN(int x, int y, int z, T omega);
  void addTemperatureBoundaryCornerPPP(int x, int y, int z, T omega);

  //  determines whether it is a corner, edge or plane boundary
  void addTemperatureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure,
                              int material, int x0, int x1, int y0, int y1,
                              int z0, int z1, T omega);
  void addTemperatureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure,
                              int material, T omega);

  void addDiffuseReflectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1, T omega, T zeta);
  void addDiffuseReflectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega, T zeta);

  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure,
                             int material, int x0, int x1, int y0, int y1, int z0, int z1);
  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material);

  void addZeroDistributionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure,
                                   int material, int x0, int x1, int y0, int y1, int z0, int z1);
  void addZeroDistributionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material);

  void addExtFieldBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure,
                           int material, int offset, int x0, int x1, int y0,
                           int y1, int z0, int z1);
  void addExtFieldBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure,
                           int material, int offset);

  virtual BlockLatticeStructure3D<T,Lattice>& getBlock();
  virtual BlockLatticeStructure3D<T,Lattice> const& getBlock() const;
private:
  template<int direction, int orientation>
  void addTemperatureBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int plane, int normal1, int normal2>
  void addTemperatureBoundaryEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int normalX, int normalY, int normalZ>
  void addTemperatureBoundaryCorner(int x, int y, int z, T omega);

private:
  BlockLatticeStructure3D<T,Lattice>& block;
  std::vector<Momenta<T,Lattice>*>  momentaVector;
  std::vector<Dynamics<T,Lattice>*> dynamicsVector;
};



} // namespace openlb


#endif
