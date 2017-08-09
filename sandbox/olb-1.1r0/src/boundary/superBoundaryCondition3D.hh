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
 * A helper for initialising 3D boundaries -- generic implementation.
 */

#ifndef SUPER_BOUNDARY_CONDITION_3D_HH
#define SUPER_BOUNDARY_CONDITION_3D_HH

#include <vector>
#include "boundaryCondition3D.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "superBoundaryCondition3D.h"
#include "core/superLattice3D.h"

namespace olb {

///////// class superBoundaryCondition3D ///////////////////////////////

template<typename T, template<typename U> class Lattice>
sOnLatticeBoundaryCondition3D<T, Lattice>::sOnLatticeBoundaryCondition3D(
  SuperLattice3D<T, Lattice>& sLattice) :
  clout(std::cout,"sOnLatticeBoundaryCondition3D"),
  _sLattice(sLattice),
  _output(false)
{
}

template<typename T, template<typename U> class Lattice>
sOnLatticeBoundaryCondition3D<T, Lattice>::sOnLatticeBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice> const& rhs) :
  clout(std::cout,"sOnLatticeBoundaryCondition3D"),
  _sLattice(rhs._sLattice),
  _output(false)
{

  _blockBCs = rhs._blockBCs;
  _overlap = rhs._overlap;
}

template<typename T, template<typename U> class Lattice>
sOnLatticeBoundaryCondition3D<T, Lattice> sOnLatticeBoundaryCondition3D<T,
                              Lattice>::operator=(sOnLatticeBoundaryCondition3D<T, Lattice> rhs)
{

  sOnLatticeBoundaryCondition3D<T, Lattice> tmp(rhs);
  return tmp;
}

template<typename T, template<typename U> class Lattice>
sOnLatticeBoundaryCondition3D<T, Lattice>::~sOnLatticeBoundaryCondition3D()
{

  for (unsigned iC = 0; iC < _blockBCs.size(); iC++) {
    delete _blockBCs[iC];
  }
}


template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addVelocityBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega)
{

  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->addVelocityBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, omega);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addSlipBoundary(
  SuperGeometry3D<T>& superGeometry, int material)
{

  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->addSlipBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addTemperatureBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega)
{

  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _ADblockBCs[iCloc]->addTemperatureBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, omega);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addDiffuseReflectionBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega, T zeta)
{

  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _ADblockBCs[iCloc]->addDiffuseReflectionBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, omega, zeta);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addConvectionBoundary(
  SuperGeometry3D<T>& superGeometry, int material)
{
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _ADblockBCs[iCloc]->addConvectionBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addZeroDistributionBoundary(
  SuperGeometry3D<T>& superGeometry, int material)
{
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _ADblockBCs[iCloc]->addZeroDistributionBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addExtFieldBoundary(
  SuperGeometry3D<T>& superGeometry, int material, int offset)
{
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _ADblockBCs[iCloc]->addExtFieldBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, offset);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addPressureBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega)
{

  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->addPressureBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, omega);
  }
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addConvectionBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega, T* uAv)
{

  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->addConvectionBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, omega, uAv);
  }
  addPoints2CommBC(superGeometry, material);
}


template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material)
{

  if (_overlap != 0) {
    int nC = _sLattice.getLoadBalancer().size();
    for (int iCloc = 0; iCloc < nC; iCloc++) {
      int nX = superGeometry.getBlockGeometry(iCloc).getNx();
      int nY = superGeometry.getBlockGeometry(iCloc).getNy();
      int nZ = superGeometry.getBlockGeometry(iCloc).getNz();

      for (int iX=-_overlap; iX<nX+_overlap; iX++) {
        for (int iY=-_overlap; iY<nY+_overlap; iY++) {
          for (int iZ=-_overlap; iZ<nZ+_overlap; iZ++) {
            if (iX < 0 || iX > nX - 1 ||
                iY < 0 || iY > nY - 1 ||
                iZ < 0 || iZ > nZ - 1 ) {
              int found = false;
              if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY,iZ)!=0) {
                for (int iXo=-_overlap; iXo<=_overlap; iXo++) {
                  for (int iYo=-_overlap; iYo<=_overlap; iYo++) {
                    for (int iZo=-_overlap; iZo<=_overlap; iZo++) {
                      int nextX = iXo + iX;
                      int nextY = iYo + iY;
                      int nextZ = iZo + iZ;
                      if (superGeometry.getBlockGeometry(iCloc).getMaterial(nextX,nextY,nextZ)==material && nextX>=-_overlap && nextY>=-_overlap && nextZ>=-_overlap && nextX<nX+_overlap && nextY<nY+_overlap && nextZ<nZ+_overlap) {
                        _sLattice.get_commBC().add_cell(iCloc, iX, iY, iZ);
                        //std::cout << "found:" <<iX<<"/"<<iY<<"/"<<iZ<<"::"<<nX<<"/"<<nY<<"/"<<nZ<<std::endl;
                        found = true;
                      }
                      if (found) {
                        break;
                      }
                    }
                    if (found) {
                      break;
                    }
                  }
                  if (found) {
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

////////////////// Factory functions //////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createLocalBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(0);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition3D<T, Lattice>* blockBC =
      createLocalBoundaryCondition3D(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createInterpBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition3D<T, Lattice>* blockBC =
      createInterpBoundaryCondition3D<T,Lattice,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createExtFdBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, Lattice>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition3D<T, Lattice>* blockBC =
      createExtendedFdBoundaryCondition3D<T,Lattice,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

//////////////// Output functions //////////////////////////////////
template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::outputOn()
{
  _output = true;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOn();
  }
}

template<typename T, template<typename U> class Lattice>
void sOnLatticeBoundaryCondition3D<T, Lattice>::outputOff()
{
  _output = false;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOff();
  }
}

} // namespace olb

#endif
