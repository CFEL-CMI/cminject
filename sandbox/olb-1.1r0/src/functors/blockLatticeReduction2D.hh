/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Albert Mink, Mathias J. Krause
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
 * The description of a blockLatticeReduction -- generic implementation.
 */

#ifndef BLOCK_LATTICE_REDUCTION_2D_HH
#define BLOCK_LATTICE_REDUCTION_2D_HH

#include<limits>
#include<cmath>

#include "functors/blockLatticeReduction2D.h"
#include "utilities/vectorHelpers.h"
#include "functors/interpolationF2D.h"
#include "functors/interpolationF3D.h"
#include "communication/mpiManager.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

// reduction
template <typename T, template <typename U> class DESCRIPTOR >
BlockLatticeReduction2D<T, DESCRIPTOR>::BlockLatticeReduction2D
(SuperLatticeF2D<T, DESCRIPTOR>& f, int resolution)
  : BlockDataF2D<T,T>( 1,1,1 ), _f(f), _origin(f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getOrigin()), _resolution(resolution),
    _h( _f.getSuperLattice().getCuboidGeometry().getMinDeltaR() ),
    _nx(_f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getNx()),
    _ny(_f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getNy()),
    clout(std::cout, "BlockLatticeReduction")
{
  this->getName() = "planeReduction(" + _f.getName() + ")";

  // [!] check whether dim(Functor) != 1
  if ( _f.getTargetDim() != 1 ) {
    clout << "Error: Functor targetDim is not 1. " << std::endl;
    exit(-1);
  }
  _origin[0] -= 10*std::numeric_limits<T>::epsilon();
  _origin[1] -= 10*std::numeric_limits<T>::epsilon();

  // changes _h, _nx, _ny if needed
  updateToWantedResolution();
  // constructs block data fields
  _tmpBlockData = new BlockData2D<T,T>( _nx, _ny );
  this->_blockData = *_tmpBlockData;
  // first update of data
  update();
}


template <typename T, template <typename U> class DESCRIPTOR >
BlockLatticeReduction2D<T, DESCRIPTOR>::~BlockLatticeReduction2D()
{
  delete _tmpBlockData;
}


template <typename T, template <typename U> class DESCRIPTOR >
void BlockLatticeReduction2D<T, DESCRIPTOR>::update()
{
  // provides physical operator()
  _f.getSuperLattice().communicate();
  AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> analyticalF( _f );

#ifdef PARALLEL_MODE_MPI
  BlockData2D<T,T> tmpBlockData( _nx, _ny );
#endif

  // -_iVoxelNmb to +_iVoxelNmb to ensure that all possible values are captured
  // since origin can be an arbitrary point in geometry
  for ( int iX = 0; iX < _nx; ++iX ) {
    for ( int iY = 0; iY < _ny; ++iY ) {
      // [!] physical units
      std::vector<T> vTmp( _origin);
      vTmp[0] += double(iX)*_h;
      vTmp[1] += double(iY)*_h;
      // store default value
      this->_blockData.get(iX, iY, 0) = T();
      // parallelization, get Cuboid Nmb out of coordinate
      int iC = 0;
      if ( _f.getSuperLattice().getCuboidGeometry().getC(vTmp, iC)  &&
           _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().checkInters(vTmp[0],vTmp[0],vTmp[1],vTmp[1]) ) {
        int rankiC = _f.getSuperLattice().getLoadBalancer().rank(iC);
        T tmp[_f.getTargetDim()];
        T vTmp2[vTmp.size()];
        for ( unsigned i = 0; i < vTmp.size(); ++i) {
          vTmp2[i] = vTmp[i];
        }
        if (analyticalF(tmp,vTmp2)) {
          if ( singleton::mpi().getRank() == rankiC ) {
            // store functor value
#ifdef PARALLEL_MODE_MPI
            tmpBlockData.get( iX, iY, 0) += tmp[0];
#else
            this->_blockData.get( iX, iY, 0) += tmp[0];
#endif
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduce(tmpBlockData, this->getBlockData(), MPI_SUM);
#endif
}

template <typename T, template <typename U> class DESCRIPTOR >
BlockStructure2D& BlockLatticeReduction2D<T, DESCRIPTOR>::getBlockStructure()
{
  return this->_blockData;
}

template <typename T, template <typename U> class DESCRIPTOR >
void BlockLatticeReduction2D<T, DESCRIPTOR>::updateToWantedResolution()
{
  if (_resolution>0) {
    if (_nx > _ny) {
      T newH = _nx*_h/(T)_resolution;
      _nx = _resolution;
      _ny = (int)(_ny*_h/newH) + 1;
      _h = newH;
    } else {
      T newH = _ny*_h/(T)_resolution;
      _ny = _resolution;
      _nx = (int)(_nx*_h/newH) + 1;
      _h = newH;
    }
  }
}

//////////////////// 3D Version ////////////////////////

// reduction
template <typename T, template <typename U> class DESCRIPTOR >
BlockLatticeReduction3D<T, DESCRIPTOR>::BlockLatticeReduction3D
(SuperLatticeF3D<T, DESCRIPTOR>& f, Vector<T,3>& u, Vector<T,3>& v,
 int resolution, T const origin[3])
  : BlockDataF2D<T,T>( 1,1,1 ), _f(f), _resolution(resolution),
    _h( _f.getSuperLattice().getCuboidGeometry().getMinDeltaR() ),
    clout(std::cout, "BlockLatticeReduction")
{
  this->getName() = "planeReduction(" + _f.getName() + ")";

  // [!] check whether dim(Functor) != 1
  if ( _f.getTargetDim() != 1 ) {
    clout << "Error: Functor targetDim is not 1. " << std::endl;
    exit(-1);
  }
  // sets default origin to the center of the guboid geometry
  constructOrigin(origin);
  // scale vectors which span the plane to be of lenght _h in Si units
  _u = u;
  _u.normalize(_h);
  _v = v;
  _v.normalize(_h);
  // computes max possible distance (sets maxLatticeDistance)
  int maxLatticeDistance = computeMaxLatticeDistance();
  // computes _origin and _nx, _ny such that the cuboid is right inside cuboid geomety and not too big (sets _origin, _nx, _ny)
  constructCuboid(maxLatticeDistance);
  // changes _h, _nx, _ny if needed
  updateToWantedResolution();
  // constructs block data fields
  _tmpBlockData = new BlockData2D<T,T>( _nx, _ny );
  this->_blockData = *_tmpBlockData;
  // first update of data
  update();
}
/*
template <typename T, template <typename U> class DESCRIPTOR >
BlockLatticeReduction3D<T, DESCRIPTOR>::BlockLatticeReduction3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, T const normal[3], int resolution, T const origin[3])
  : BlockF2D<T>( *_blockData, 1 ), _f(f), _origin(3,T()),
    _resolution(resolution),
    _h( _f.getSuperLattice().getCuboidGeometry().getMinDeltaR() ), _u(3,T()),
    _v(3,T()), clout(std::cout, "BlockLatticeReduction")
{
  this->getName() = "planeReduction(" + _f.getName() + ")";

  // [!] check whether dim(Functor) != 1
  if ( _f.getTargetDim() != 1 ) {
    clout << "Error: Functor targetDim is not 1. " << std::endl;
    exit(-1);
  }

  // sets default origin to the center of the guboid geometry
  constructOrigin(origin);

  // scale vectors which span the plane to be of lenght _h in Si units
  if (!(normal[0]*normal[1]*normal[2]) ) {
    if (!(normal[0]) ) {
      _u[0] = T(_h);
      _u[1] = T();
      _u[2] = T();
    } else if (!(normal[1]) ) {
      _u[1] = T(_h);
      _u[0] = T();
      _u[2] = T();
    } else if (!(normal[2]) ) {
      _u[2] = T(_h);
      _u[0] = T();
      _u[1] = T();
    }
  } else {
    _u[0] = T();
    _u[1] = normal[1];
    _u[2] = -normal[2];
    _u = normalize(_u)*_h;
  }
  _v[0] = normal[0];
  _v[1] = normal[1];
  _v[2] = normal[2];
  _v = normalize(crossProduct3D(_u,_v))*_h;

  // computes max possible distance (sets maxLatticeDistance)
  int maxLatticeDistance = computeMaxLatticeDistance();

  // computes _origin and _nx, _ny such that the cuboid is right inside cuboid geomety and not too big (sets _origin, _nx, _ny)
  constructCuboid(maxLatticeDistance);

  // changes _h, _nx, _ny if needed
  updateToWantedResolution();

  // constructs block data fields
  _blockData = new BlockData2D<T>( _nx, _ny );
  _tmpBlockData = new BlockData2D<T>( _nx, _ny );

  // first update of data
  update();
}
*/
template <typename T, template <typename U> class DESCRIPTOR >
BlockLatticeReduction3D<T, DESCRIPTOR>::BlockLatticeReduction3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, T const normalX, T const normalY, T const normalZ, int resolution)
  : BlockDataF2D<T,T>( 1,1,1 ), _f(f), _resolution(resolution),
    _h( _f.getSuperLattice().getCuboidGeometry().getMinDeltaR() ),
    clout(std::cout, "BlockLatticeReduction")
{
  this->getName() = "planeReduction(" + _f.getName() + ")";

  // [!] check whether dim(Functor) != 1
  if ( _f.getTargetDim() != 1 ) {
    clout << "Error: Functor targetDim is not 1. " << std::endl;
    exit(-1);
  }

  // sets default origin to the center of the guboid geometry
  constructOrigin();

  // scale vectors which span the plane to be of lenght _h in Si units
  Vector<T,3> normal = {normalX, normalY, normalZ};
  if ( util::nearZero(normal[0]*normal[1]*normal[2]) ) {
    if ( util::nearZero(normal[0]) ) {
      _u = {T(_h), T(), T()};
    } else if ( util::nearZero(normal[1]) ) {
      _u = {T(), T(_h), T()};
    } else if ( util::nearZero(normal[2]) ) {
      _u = {T(), T(), T(_h)};
    }
  } else {
    _u = {T(), normal[1], -normal[2]};
    _u.normalize(_h);
  }
  _v = normal;
  _v = crossProduct3D(_u,_v);
  _v.normalize(_h);

  // computes max possible distance (sets maxLatticeDistance)
  int maxLatticeDistance = computeMaxLatticeDistance();
  // computes _origin and _nx, _ny such that the cuboid is right inside cuboid geomety and not too big (sets _origin, _nx, _ny)
  constructCuboid(maxLatticeDistance);
  // changes _h, _nx, _ny if needed
  updateToWantedResolution();
  // constructs block data fields
  _tmpBlockData = new BlockData2D<T,T>( _nx, _ny );
  this->_blockData = *_tmpBlockData;
  // first update of data
  update();
}

template <typename T, template <typename U> class DESCRIPTOR >
BlockLatticeReduction3D<T, DESCRIPTOR>::BlockLatticeReduction3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, T const normalX, T const normalY, T const normalZ,
  T const originX, T const originY, T const originZ, int resolution)
  : BlockDataF2D<T,T>( 1,1,1 ), _f(f), _resolution(resolution),
    _h( _f.getSuperLattice().getCuboidGeometry().getMinDeltaR() ),
    clout(std::cout, "BlockLatticeReduction")
{
  this->getName() = "planeReduction(" + _f.getName() + ")";

  // [!] check whether dim(Functor) != 1
  if ( _f.getTargetDim() != 1 ) {
    clout << "Error: Functor targetDim is not 1. " << std::endl;
    exit(-1);
  }

  // sets default origin to the center of the guboid geometry
  T origin[3] = {originX, originY, originZ};
  constructOrigin(origin);

  // scale vectors which span the plane to be of lenght _h in Si units
  Vector<T,3> normal = {normalX, normalY, normalZ};
  if ( util::nearZero(normal[0]*normal[1]*normal[2]) ) {
    if ( util::nearZero(normal[0]*normal[1]*normal[2]) ) {
      if ( util::nearZero(normal[0]) ) {
        _u = {T(_h), T(), T()};
      } else if ( util::nearZero(normal[1]) ) {
        _u = {T(), T(_h), T()};
      } else if ( util::nearZero(normal[2]) ) {
        _u = {T(), T(), T(_h)};
      }
    } else {
      _u = {T(), normal[1], -normal[2]};
      _u.normalize(_h);
    }
    _v = normal;
    _v = crossProduct3D(_u,_v);
    _v.normalize(_h);
    // computes max possible distance (sets maxLatticeDistance)
    int maxLatticeDistance = computeMaxLatticeDistance();

    // computes _origin and _nx, _ny such that the cuboid is right inside cuboid geomety and not too big (sets _origin, _nx, _ny)
    constructCuboid(maxLatticeDistance);
    // changes _h, _nx, _ny if needed
    updateToWantedResolution();
    // constructs block data fields
    _tmpBlockData = new BlockData2D<T,T>( _nx, _ny );
    this->_blockData = *_tmpBlockData;
    // first update of data
    update();
  }
}

template <typename T, template <typename U> class DESCRIPTOR >
BlockLatticeReduction3D<T, DESCRIPTOR>::~BlockLatticeReduction3D()
{
  delete _tmpBlockData;
}


template <typename T, template <typename U> class DESCRIPTOR >
void BlockLatticeReduction3D<T, DESCRIPTOR>::constructOrigin(T const origin[3])
{
  if (origin == nullptr) {
    Vector<T,3> originCuboid( _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getOrigin() );
    Vector<int,3> extend( _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getExtend() );
    T deltaR = _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getDeltaR();
    _origin[0] = (originCuboid[0] + 0.5 * deltaR * extend[0]);
    _origin[1] = (originCuboid[1] + 0.5 * deltaR * extend[1]);
    _origin[2] = (originCuboid[2] + 0.5 * deltaR * extend[2]);
    _origin[0] -= 2*std::numeric_limits<T>::epsilon()*fabs(_origin[0]);
    _origin[1] -= 2*std::numeric_limits<T>::epsilon()*fabs(_origin[1]);
    _origin[2] -= 2*std::numeric_limits<T>::epsilon()*fabs(_origin[2]);
  } else {
    _origin[0] = origin[0] - 2*std::numeric_limits<T>::epsilon()*fabs(origin[0]);
    _origin[1] = origin[1] - 2*std::numeric_limits<T>::epsilon()*fabs(origin[1]);
    _origin[2] = origin[2] - 2*std::numeric_limits<T>::epsilon()*fabs(origin[2]);
  }
}

template <typename T, template <typename U> class DESCRIPTOR >
int BlockLatticeReduction3D<T, DESCRIPTOR>::computeMaxLatticeDistance()
{
  Vector<T,3> originCuboid( _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getOrigin() );
  Vector<int,3> extend( _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getExtend() );
  T deltaR = _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().getDeltaR();
  T maxPhysDistance = T();
  T tmp;
  for (int iDim = 0; iDim < 3; ++iDim) {
    tmp = (originCuboid[iDim] - _origin[iDim])*(originCuboid[iDim] - _origin[iDim]);
    if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
    }
    tmp = (originCuboid[iDim] + extend[iDim]*deltaR - _origin[iDim])*(originCuboid[iDim] + extend[iDim]*deltaR - _origin[iDim]);
    if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
    }
  }
  return sqrt(maxPhysDistance)/_h + 1;
}


template <typename T, template <typename U> class DESCRIPTOR >
void BlockLatticeReduction3D<T, DESCRIPTOR>::constructCuboid(int maxLatticeDistance)
{
  int iC;
  int minX = -maxLatticeDistance;
  int maxX = maxLatticeDistance;
  bool found = false;
  for ( int iX = -maxLatticeDistance; iX < maxLatticeDistance; ++iX ) {
    for ( int iY = -maxLatticeDistance; iY < maxLatticeDistance; ++iY ) {
      std::vector<T> physR(3,T());
      physR[0] = _origin[0] + double(iX)*_u[0] + double(iY)*_v[0];
      physR[1] = _origin[1] + double(iX)*_u[1] + double(iY)*_v[1];
      physR[2] = _origin[2] + double(iX)*_u[2] + double(iY)*_v[2];

      if ( _f.getSuperLattice().getCuboidGeometry().getC(physR, iC) ) {
        minX = iX;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }
  found = false;
  for ( int iX = maxLatticeDistance; iX > -maxLatticeDistance; --iX ) {
    for ( int iY = -maxLatticeDistance; iY < maxLatticeDistance; ++iY ) {
      std::vector<T> physR(3,T());
      physR[0] = _origin[0] + double(iX)*_u[0] + double(iY)*_v[0];
      physR[1] = _origin[1] + double(iX)*_u[1] + double(iY)*_v[1];
      physR[2] = _origin[2] + double(iX)*_u[2] + double(iY)*_v[2];

      if ( _f.getSuperLattice().getCuboidGeometry().getC(physR, iC) ) {
        maxX = iX;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }

  int minY = -maxLatticeDistance;
  int maxY = maxLatticeDistance;
  found = false;
  for ( int iY = -maxLatticeDistance; iY < maxLatticeDistance; ++iY ) {
    for ( int iX = -maxLatticeDistance; iX < maxLatticeDistance; ++iX ) {
      std::vector<T> physR(3,T());
      physR[0] = _origin[0] + double(iX)*_u[0] + double(iY)*_v[0];
      physR[1] = _origin[1] + double(iX)*_u[1] + double(iY)*_v[1];
      physR[2] = _origin[2] + double(iX)*_u[2] + double(iY)*_v[2];

      if ( _f.getSuperLattice().getCuboidGeometry().getC(physR, iC) ) {
        minY = iY;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }
  found = false;
  for ( int iY = maxLatticeDistance; iY > -maxLatticeDistance; --iY ) {
    for ( int iX = -maxLatticeDistance; iX < maxLatticeDistance; ++iX ) {
      std::vector<T> physR(3,T());
      physR[0] = _origin[0] + double(iX)*_u[0] + double(iY)*_v[0];
      physR[1] = _origin[1] + double(iX)*_u[1] + double(iY)*_v[1];
      physR[2] = _origin[2] + double(iX)*_u[2] + double(iY)*_v[2];

      if ( _f.getSuperLattice().getCuboidGeometry().getC(physR, iC) ) {
        maxY = iY;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }

  _nx = maxX - minX + 1;
  _ny = maxY - minY + 1;
  _origin = _origin + double(minX)*_u + double(minY)*_v;
}

template <typename T, template <typename U> class DESCRIPTOR >
void BlockLatticeReduction3D<T, DESCRIPTOR>::update()
{
  // provides physical operator()
  _f.getSuperLattice().communicate();
  AnalyticalFfromSuperF3D<T> analyticalF( _f );

#ifdef PARALLEL_MODE_MPI
  BlockData2D<T,T> tmpBlockData( _nx, _ny );
#endif

  // -_iVoxelNmb to +_iVoxelNmb to ensure that all possible values are captured
  // since origin can be an arbitrary point in geometry
  for ( int iX = 0; iX < _nx; ++iX ) {
    for ( int iY = 0; iY < _ny; ++iY ) {
      // [!] physical units
      std::vector<T> vTmp(3,T());
      vTmp[0] = _origin[0] + double(iX)*_u[0] + double(iY)*_v[0];
      vTmp[1] = _origin[1] + double(iX)*_u[1] + double(iY)*_v[1];
      vTmp[2] = _origin[2] + double(iX)*_u[2] + double(iY)*_v[2];

      // store default value
      this->_blockData.get(iX, iY, 0) = T();
      // parallelization, get Cuboid Nmb out of coordinate
      int iC;
      if ( _f.getSuperLattice().getCuboidGeometry().getC(vTmp, iC) &&
           _f.getSuperLattice().getCuboidGeometry().getMotherCuboid().checkInters(vTmp[0],vTmp[0],vTmp[1],vTmp[1],vTmp[2],vTmp[2]) ) {
        int rankiC = _f.getSuperLattice().getLoadBalancer().rank(iC);
        T tmp[_f.getTargetDim()];
        T vTmp2[vTmp.size()];
        for (unsigned i = 0; i < vTmp.size(); ++i) {
          vTmp2[i] = vTmp[i];
        }
        if (analyticalF(tmp,vTmp2)) {
          if ( singleton::mpi().getRank() == rankiC ) {
            // store functor value
#ifdef PARALLEL_MODE_MPI
            tmpBlockData.get( iX, iY, 0) += tmp[0];
#else
            this->_blockData.get( iX, iY, 0) += tmp[0];
#endif
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduce(tmpBlockData, this->getBlockData(), MPI_SUM);
#endif
}


template <typename T, template <typename U> class DESCRIPTOR >
BlockStructure2D& BlockLatticeReduction3D<T, DESCRIPTOR>::getBlockStructure()
{
  return this->_blockData;
}

template <typename T, template <typename U> class DESCRIPTOR >
void BlockLatticeReduction3D<T, DESCRIPTOR>::updateToWantedResolution()
{
  if (_resolution > 0) {
    if (_nx > _ny) {
      T newH = _nx*_h/(T)_resolution;
      _nx = _resolution;
      _ny = (int)(_ny*_h/newH) + 1;
      _h = newH;
    } else {
      T newH = _ny*_h/(T)_resolution;
      _ny = _resolution;
      _nx = (int)(_nx*_h/newH) + 1;
      _h = newH;
    }
  }
  _u.normalize();
  _u*=_h;
  _v.normalize();
  _v*=_h;
}

} // end namespace olb

#endif
