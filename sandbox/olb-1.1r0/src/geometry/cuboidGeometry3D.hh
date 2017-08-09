/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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
 * The description of a vector of 3D cuboid  -- generic implementation.
 */


#ifndef CUBOID_GEOMETRY_3D_HH
#define CUBOID_GEOMETRY_3D_HH


#include <iostream>
#include <math.h>
#include <algorithm>
#include <set>
#include <limits>
#include "geometry/cuboidGeometry3D.h"
#include "functors/indicator/indicatorF3D.h"
#include "communication/loadBalancer.h"


namespace olb {

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D()
  : _motherCuboid(0,0,0,0,0,0,0), _periodicityOn(false), clout(std::cout, "CuboidGeometry3D")
{
  add(_motherCuboid);
  split(0, 1);
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(T originX, T originY, T originZ, T deltaR,
                                      int nX, int nY, int nZ, int nC)
  : _motherCuboid(originX, originY, originZ, deltaR, nX, nY, nZ),
    _periodicityOn(false), clout(std::cout, "CuboidGeometry3D")
{
  add(_motherCuboid);
  split(0, nC);
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(std::vector<T> origin, T deltaR,
                                      std::vector<int> extent, int nC)
  : CuboidGeometry3D(origin[0], origin[1], origin[2], deltaR,
                     extent[0], extent[1], extent[2], nC)
{
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(IndicatorF3D<T>& indicatorF, T voxelSize, int nC)
  : _motherCuboid( indicatorF.getMin()[0],  indicatorF.getMin()[1], indicatorF.getMin()[2], voxelSize,
                   (int)((indicatorF.getMax()[0] - indicatorF.getMin()[0]) / voxelSize + 1.5),
                   (int)((indicatorF.getMax()[1] - indicatorF.getMin()[1]) / voxelSize + 1.5),
                   (int)((indicatorF.getMax()[2] - indicatorF.getMin()[2]) / voxelSize + 1.5)),
    _periodicityOn(false), clout(std::cout, "CuboidGeometry3D")
{
  add(_motherCuboid);
  split(0, nC);
  shrink(indicatorF);
}

template<typename T>
CuboidGeometry3D<T>::~CuboidGeometry3D() {};


template<typename T>
Cuboid3D<T>& CuboidGeometry3D<T>::get(int iC)
{
  return _cuboids[iC];
}

template<typename T>
Cuboid3D<T> const& CuboidGeometry3D<T>::get(int iC) const
{
  return _cuboids[iC];
}

template<typename T>
Cuboid3D<T> CuboidGeometry3D<T>::getMotherCuboid()
{
  return _motherCuboid;
}

template<typename T>
Cuboid3D<T> const& CuboidGeometry3D<T>::getMotherCuboid() const
{
  return _motherCuboid;
}

template<typename T>
void CuboidGeometry3D<T>::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ)
{
  _periodicityOn[0] = periodicityX;
  _periodicityOn[1] = periodicityY;
  _periodicityOn[2] = periodicityZ;
}


template<typename T>
int CuboidGeometry3D<T>::get_iC(T x, T y, T z, int offset) const
{
  unsigned i;
  for (i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(x, y, z, offset)) {
      return (int)i;
    }
  }
  return (int)i;
}


template<typename T>
int CuboidGeometry3D<T>::get_iC(Vector<T,3> coords, int offset) const
{
  return get_iC(coords[0], coords[1], coords[2], offset);
}

template<typename T>
int CuboidGeometry3D<T>::get_iC(T x, T y, T z, int orientationX, int orientationY,
                                int orientationZ) const
{
  unsigned i;
  for (i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(x, y, z) &&
        _cuboids[i].checkPoint(x + orientationX / _cuboids[i].getDeltaR(),
                               y + orientationY / _cuboids[i].getDeltaR(),
                               z + orientationZ / _cuboids[i].getDeltaR())) {
      return (int)i;
    }
  }
  return (int)i;
}

template<typename T>
bool CuboidGeometry3D<T>::getC(std::vector<T> physR, int& iC) const
{
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp < getNc()) {
    iC = iCtmp;
    return true;
  } else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry3D<T>::getLatticeR(int latticeR[4], const T physR[3]) const
{
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    latticeR[3] = (int)floor( (physR[2] - _cuboids[latticeR[0]].getOrigin()[2] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    return true;
  } else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry3D<T>::getFloorLatticeR(const std::vector<T>& physR, std::vector<int>& latticeR) const
{
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[3] = (int)floor( (physR[2] - _cuboids[latticeR[0]].getOrigin()[2] ) / _cuboids[latticeR[0]].getDeltaR() );
    return true;
  } else {
    return false;
  }
}


template<typename T>
void CuboidGeometry3D<T>::getPhysR(T physR[3], const int& iCglob, const int& iX, const int& iY, const int& iZ) const
{
  _cuboids[iCglob].getPhysR(physR, iX, iY, iZ);
  for (int iDim = 0; iDim < 3; iDim++) {
    if (_periodicityOn[iDim]) {
      //std::cout << iDim << _periodicityOn[iDim] <<":"<< _motherCuboid.getDeltaR()*(_motherCuboid.getExtend()[iDim]) << std::endl;
      physR[iDim] = remainder( physR[iDim] - _motherCuboid.getOrigin()[iDim]
                               + _motherCuboid.getDeltaR() * (_motherCuboid.getExtend()[iDim]) ,
                               _motherCuboid.getDeltaR() * (_motherCuboid.getExtend()[iDim]));
      // solving the rounding error problem for double
      if ( physR[iDim]*physR[iDim] < 0.001 * _motherCuboid.getDeltaR()*_motherCuboid.getDeltaR() ) {
        if ( physR[iDim] > 0 ) {
          physR[iDim] = _motherCuboid.getDeltaR() * (_motherCuboid.getExtend()[iDim]);
        } else {
          physR[iDim] = T();
        }
      }
      // make it to mod instead remainer
      if ( physR[iDim] < 0 ) {
        physR[iDim] += _motherCuboid.getDeltaR() *( _motherCuboid.getExtend()[iDim]);
      }
      // add origin
      physR[iDim] += _motherCuboid.getOrigin()[iDim];
    }
  }
  return;
}

template<typename T>
void CuboidGeometry3D<T>::getPhysR(T physR[3], const int latticeR[4]) const
{
  getPhysR(physR, latticeR[0],  latticeR[1],  latticeR[2],  latticeR[3]);
}


template<typename T>
void CuboidGeometry3D<T>::getNeighbourhood(int cuboid, std::vector<int>& neighbours, int overlap)
{
  neighbours.clear();

  std::set<int> dummy;

  for (int iC = 0; iC < getNc(); iC++) {
    if (cuboid == iC) {
      continue;
    }
    T globX = get(iC).getOrigin()[0];
    T globY = get(iC).getOrigin()[1];
    T globZ = get(iC).getOrigin()[2];
    T nX = get(iC).getNx();
    T nY = get(iC).getNy();
    T nZ = get(iC).getNz();
    T deltaR = get(iC).getDeltaR();
    if (get(cuboid).checkInters(globX - overlap * deltaR,
                                globX + (nX + overlap - 1)*deltaR,
                                globY - overlap * deltaR,
                                globY + (nY + overlap - 1)*deltaR,
                                globZ - overlap * deltaR,
                                globZ + (nZ + overlap - 1)*deltaR, overlap)) {
      //neighbours.push_back(iC);
      dummy.insert(iC);
    }

    if (_periodicityOn[0]) {
      if (get(cuboid).getOrigin()[0] + (get(cuboid).getNx() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[0]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0]-getMaxPhysR()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
      if (get(cuboid).getOrigin()[0] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[0]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0]+getMaxPhysR()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
    }

    if (_periodicityOn[1]) {
      if (get(cuboid).getOrigin()[1] + (get(cuboid).getNy() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[1]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1]-getMaxPhysR()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
      if (get(cuboid).getOrigin()[1] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[1]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1]+getMaxPhysR()[1],
                        get(cuboid).getOrigin()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
    }

    if (_periodicityOn[2]) {
      if (get(cuboid).getOrigin()[2] + (get(cuboid).getNz() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[2]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2]-getMaxPhysR()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
      if (get(cuboid).getOrigin()[2] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[2]) {
        Cuboid3D<T> cub(get(cuboid).getOrigin()[0],
                        get(cuboid).getOrigin()[1],
                        get(cuboid).getOrigin()[2]+getMaxPhysR()[2],
                        get(cuboid).getDeltaR(),
                        get(cuboid).getNx(),
                        get(cuboid).getNy(),
                        get(cuboid).getNz());
        if (cub.checkInters(globX - overlap * deltaR,
                            globX + (nX + overlap - 1)*deltaR,
                            globY - overlap * deltaR,
                            globY + (nY + overlap - 1)*deltaR,
                            globZ - overlap * deltaR,
                            globZ + (nZ + overlap - 1)*deltaR, overlap)) {
          dummy.insert(iC);
        }
      }
    }

  }
  std::set<int>::iterator it = dummy.begin();
  for (; it != dummy.end(); ++it) {
    neighbours.push_back(*it);
  }
}

template<typename T>
int CuboidGeometry3D<T>::getNc() const
{
  return _cuboids.size();
}

template<typename T>
T CuboidGeometry3D<T>::getMinRatio() const
{
  T minRatio = 1.;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if ((T)_cuboids[i].getNx() / (T)_cuboids[i].getNy() < minRatio) {
      minRatio = (T)_cuboids[i].getNx() / (T)_cuboids[i].getNy();
    }
    if ((T)_cuboids[i].getNy() / (T)_cuboids[i].getNz() < minRatio) {
      minRatio = (T)_cuboids[i].getNy() / (T)_cuboids[i].getNz();
    }
    if ((T)_cuboids[i].getNz() / (T)_cuboids[i].getNx() < minRatio) {
      minRatio = (T)_cuboids[i].getNz() / (T)_cuboids[i].getNx();
    }
  }
  return minRatio;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxRatio() const
{
  T maxRatio = 1.;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if ((T)_cuboids[i].getNx() / (T)_cuboids[i].getNy() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNx() / (T)_cuboids[i].getNy();
    }
    if ((T)_cuboids[i].getNy() / (T)_cuboids[i].getNz() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNy() / (T)_cuboids[i].getNz();
    }
    if ((T)_cuboids[i].getNz() / (T)_cuboids[i].getNx() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNz() / (T)_cuboids[i].getNx();
    }
  }
  return maxRatio;
}

template<typename T>
Vector<T,3> CuboidGeometry3D<T>::getMinPhysR() const
{
  Vector<T,3> output (_cuboids[0].getOrigin());
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getOrigin()[0] < output[0]) {
      output[0] = _cuboids[i].getOrigin()[0];
    }
    if (_cuboids[i].getOrigin()[1] < output[1]) {
      output[1] = _cuboids[i].getOrigin()[1];
    }
    if (_cuboids[i].getOrigin()[2] < output[2]) {
      output[2] = _cuboids[i].getOrigin()[2];
    }
  }
  return output;
}

template<typename T>
Vector<T,3> CuboidGeometry3D<T>::getMaxPhysR() const
{
  Vector<T,3> output (_cuboids[0].getOrigin());
  output[0] += _cuboids[0].getNx()*_cuboids[0].getDeltaR();
  output[1] += _cuboids[0].getNy()*_cuboids[0].getDeltaR();
  output[2] += _cuboids[0].getNz()*_cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getOrigin()[0] + _cuboids[i].getNx()*_cuboids[i].getDeltaR() > output[0]) {
      output[0] = _cuboids[i].getOrigin()[0] + _cuboids[i].getNx()*_cuboids[i].getDeltaR();
    }
    if (_cuboids[i].getOrigin()[1] + _cuboids[i].getNy()*_cuboids[i].getDeltaR() > output[1]) {
      output[1] = _cuboids[i].getOrigin()[1] + _cuboids[i].getNy()*_cuboids[i].getDeltaR();
    }
    if (_cuboids[i].getOrigin()[2] + _cuboids[i].getNz()*_cuboids[i].getDeltaR() > output[2]) {
      output[2] = _cuboids[i].getOrigin()[2] + _cuboids[i].getNz()*_cuboids[i].getDeltaR();
    }
  }
  return output;
}

template<typename T>
T CuboidGeometry3D<T>::getMinPhysVolume() const
{
  T minVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() < minVolume) {
      minVolume = _cuboids[i].getPhysVolume();
    }
  }
  return minVolume;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxPhysVolume() const
{
  T maxVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() > maxVolume) {
      maxVolume = _cuboids[i].getPhysVolume();
    }
  }
  return maxVolume;
}

template<typename T>
int CuboidGeometry3D<T>::getMinLatticeVolume() const
{
  int minNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() < minNodes) {
      minNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return minNodes;
}

template<typename T>
int CuboidGeometry3D<T>::getMaxLatticeVolume() const
{
  int maxNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() > maxNodes) {
      maxNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return maxNodes;
}

template<typename T>
T CuboidGeometry3D<T>::getMinDeltaR() const
{
  T minDelta = _cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getDeltaR() < minDelta) {
      minDelta = _cuboids[i].getDeltaR();
    }
  }
  return minDelta;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxDeltaR() const
{
  T maxDelta = _cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getDeltaR() > maxDelta) {
      maxDelta = _cuboids[i].getDeltaR();
    }
  }
  return maxDelta;
}


template<typename T>
bool CuboidGeometry3D<T>::operator==(CuboidGeometry3D<T>& rhs)
{
  return     _motherCuboid == rhs._motherCuboid
             && _periodicityOn == rhs._periodicityOn
             &&       _cuboids == rhs._cuboids;
}



template<typename T>
void CuboidGeometry3D<T>::add(Cuboid3D<T> cuboid)
{

  _cuboids.push_back(cuboid);
}

template<typename T>
void CuboidGeometry3D<T>::remove(int iC)
{

  _cuboids.erase(_cuboids.begin() + iC);
}


template<typename T>
void CuboidGeometry3D<T>::remove(IndicatorF3D<T>& indicatorF)
{

  //IndicatorIdentity3D<T> tmpIndicatorF(indicatorF);

  std::vector<bool> allZero;
  int latticeR[4];
  T physR[3];
  for (unsigned iC = 0; iC < _cuboids.size(); iC++) {
    latticeR[0] = iC;
    allZero.push_back(true);
    for (int iX = 0; iX < _cuboids[iC].getNx(); iX++) {
      for (int iY = 0; iY < _cuboids[iC].getNy(); iY++) {
        for (int iZ = 0; iZ < _cuboids[iC].getNz(); iZ++) {
          latticeR[1] = iX;
          latticeR[2] = iY;
          latticeR[3] = iZ;
          getPhysR(physR,latticeR);
          bool inside[1];
          indicatorF(inside,physR);
          if (inside[0]) {
            allZero[iC] = 0;
          }
        }
      }
    }
  }
  for (int iC = _cuboids.size() - 1; iC >= 0; iC--) {
    if (allZero[iC] ) {
      remove(iC);
    }
  }
}

template<typename T>
void CuboidGeometry3D<T>::shrink(IndicatorF3D<T>& indicatorF)
{

  //IndicatorIdentity3D<T> tmpIndicatorF(indicatorF);
  int newX, newY, newZ, maxX, maxY, maxZ;
  int nC = getNc();
  int latticeR[4];
  T physR[3];
  bool inside[1];
  for (int iC = nC - 1; iC >= 0; iC--) {
    latticeR[0] = iC;
    int fullCells = 0;
    int xN = get(iC).getNx();
    int yN = get(iC).getNy();
    int zN = get(iC).getNz();
    maxX = 0;
    maxY = 0;
    maxZ = 0;
    newX = xN - 1;
    newY = yN - 1;
    newZ = zN - 1;
    for (int iX = 0; iX < xN; iX++) {
      for (int iY = 0; iY < yN; iY++) {
        for (int iZ = 0; iZ < zN; iZ++) {
          latticeR[1] = iX;
          latticeR[2] = iY;
          latticeR[3] = iZ;
          getPhysR(physR,latticeR);
          indicatorF(inside,physR);
          if (inside[0]) {
            fullCells++;
            maxX = std::max(maxX, iX);
            maxY = std::max(maxY, iY);
            maxZ = std::max(maxZ, iZ);
            newX = std::min(newX, iX);
            newY = std::min(newY, iY);
            newZ = std::min(newZ, iZ);
          }
        }
      }
    }
    //    if (maxX+2 < xN) maxX+=2; else if (maxX+1 < xN) maxX+=1;
    //    if (maxY+2 < yN) maxY+=2; else if (maxY+1 < yN) maxY+=1;
    //    if (maxZ+2 < zN) maxZ+=2; else if (maxZ+1 < zN) maxZ+=1;
    //
    //    if (newX-2 >= 0) newX-=2; else if (newX-1 >= 0) newX-=1;
    //    if (newY-2 >= 0) newY-=2; else if (newY-1 >= 0) newY-=1;
    //    if (newZ-2 >= 0) newZ-=2; else if (newZ-1 >= 0) newZ-=1;

    if (fullCells > 0) {
      get(iC).setWeight(fullCells);
      _cuboids[iC].resize(newX, newY, newZ, maxX - newX + 1, maxY - newY + 1, maxZ - newZ + 1);
    } else {
      remove(iC);
    }
  }
  // shrink mother cuboid
  Vector<T,3> minPhysR = getMinPhysR();
  Vector<T,3> maxPhysR = getMaxPhysR();
  T minDelataR = getMinDeltaR();
  _motherCuboid = Cuboid3D<T>(minPhysR[0], minPhysR[1], minPhysR[2], minDelataR,
                              (int)((maxPhysR[0]-minPhysR[0])/minDelataR + 0.5) ,
                              (int)((maxPhysR[1]-minPhysR[1])/minDelataR + 0.5),
                              (int)((maxPhysR[2]-minPhysR[2])/minDelataR + 0.5));
}


template<typename T>
void CuboidGeometry3D<T>::split(int iC, int p)
{

  Cuboid3D<T> temp(_cuboids[iC].getOrigin()[0], _cuboids[iC].getOrigin()[1],
                   _cuboids[iC].getOrigin()[2],  _cuboids[iC].getDeltaR(),
                   _cuboids[iC].getNx(), _cuboids[iC].getNy(), _cuboids[iC].getNz());
  temp.divide(p, _cuboids);
  remove(iC);
}


template<typename T>
void CuboidGeometry3D<T>::swap(CuboidGeometry3D<T>& rhs)
{
  std::swap(this->_cuboids, rhs._cuboids);
  std::swap(this->_motherCuboid, rhs._motherCuboid);
  std::swap(this->_periodicityOn[0], rhs._periodicityOn[0]);
  std::swap(this->_periodicityOn[1], rhs._periodicityOn[1]);
  std::swap(this->_periodicityOn[2], rhs._periodicityOn[2]);
  std::swap(this->clout, rhs.clout);
}

template<typename T>
void CuboidGeometry3D<T>::swapCuboids(std::vector< Cuboid3D<T> >& cuboids)
{
  _cuboids.swap(cuboids);
}

template<typename T>
void CuboidGeometry3D<T>::replaceCuboids(std::vector< Cuboid3D<T> >& cuboids)
{
  this->_cuboids.clear();
  for ( unsigned iC = 0; iC < cuboids.size(); iC++) {
    add(cuboids[iC]);
  }
}


template<typename T>
void CuboidGeometry3D<T>::setWeights(IndicatorF3D<T>& indicatorF)
{

  //IndicatorIdentity3D<T> tmpIndicatorF(indicatorF);
  int xN, yN, zN;
  int nC = getNc();
  int latticeR[4];
  T physR[3];
  for ( int iC = nC - 1; iC >= 0; iC--) { // assemble neighbourhood information
    latticeR[0] = iC;
    xN  = get(iC).getNx();
    yN  = get(iC).getNy();
    zN  = get(iC).getNz();
    int fullCells = 0;
    for (int iX = 0; iX < xN; iX++) {
      for (int iY = 0; iY < yN; iY++) {
        for (int iZ = 0; iZ < zN; iZ++) {
          latticeR[1] = iX;
          latticeR[2] = iY;
          latticeR[3] = iZ;
          getPhysR(physR,latticeR);
          bool inside[1];
          indicatorF(inside,physR);
          if (inside[0]) {
            fullCells++;
          }
        }
      }
    }
    if (fullCells > 0) {
      get(iC).setWeight(fullCells);
    } else {
      remove(iC);
    }
  }
}


template<typename T>
size_t CuboidGeometry3D<T>::getNblock() const
{
  return   1  // _periodicityOn
           + _motherCuboid.getNblock() // _motherCuboid
           + _cuboids.size() > 0 ? 1 + _cuboids.size() * _cuboids[0].getNblock() : 0; // _cuboids;
}


template<typename T>
size_t CuboidGeometry3D<T>::getSerializableSize() const
{
  return   3 * sizeof(bool)  // _periodicityOn
           + _motherCuboid.getSerializableSize() // _motherCuboid
           + (_cuboids.size() > 0 ?
              sizeof(size_t) + _cuboids.size() * _cuboids[0].getSerializableSize() :
              0); // _cuboids;
}

template<typename T>
bool* CuboidGeometry3D<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  size_t sizeBufferIndex = 0;
  bool* dataPtr = nullptr;

  registerVar<bool>                           (iBlock, sizeBlock, currentBlock, dataPtr, _periodicityOn[0], 3);
  registerSerializableOfConstSize             (iBlock, sizeBlock, currentBlock, dataPtr, _motherCuboid, loadingMode);
  registerStdVectorOfSerializablesOfConstSize (iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr,
      _cuboids, loadingMode);

  return dataPtr;
}


template<typename T>
void CuboidGeometry3D<T>::print() const
{
  clout << "---Cuboid Stucture Statistics---" << std::endl;
  clout << " Number of Cuboids: " << "\t" << getNc() << std::endl;
  clout << " Delta (min): " << "\t" << "\t" << getMinDeltaR() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxDeltaR() << std::endl;
  clout << " Ratio (min): " << "\t" << "\t" << getMinRatio() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxRatio() << std::endl;
  clout << " Nodes (min): " << "\t" << "\t" << getMinLatticeVolume() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxLatticeVolume() << std::endl;
  clout << "--------------------------------" << std::endl;
}

template<typename T>
void CuboidGeometry3D<T>::printExtended()
{
  clout << "Mothercuboid :" << std::endl;
  getMotherCuboid().print();

  for (int iC = 0; iC < getNc(); iC++) {
    clout << "Cuboid #" << iC << ": " << std::endl;
    get(iC).print();
  }
}

template<typename T>
void CuboidGeometry3D<T>::writeToExistingFile(std::string completeFileName, LoadBalancer<T>& loadBalancer)
{
  std::ofstream fout;
  if ( singleton::mpi().isMainProcessor() ) {

    // Open File
    fout.open(completeFileName.c_str(), std::ios::app);
    if (!fout) {
      clout << "Error: could not open " << completeFileName << std::endl;
    }

    // --- Preamble --- //
    fout << "<CuboidGeometry dimension=\"3\" " << _cuboidParameters(getMotherCuboid()) << ">\n";

    // TODO: Move Cuboid XML Serialization to Cuboid3D class
    for (int iC = 0; iC < getNc(); ++iC) {
      fout << "<Cuboid " << _cuboidParameters(get(iC)) << " />\n";
    }

    fout << "</CuboidGeometry>\n";

    // Close File
    fout.close();
  }
}


template<typename T>
void CuboidGeometry3D<T>::writeToFile(std::string fileName, LoadBalancer<T>& loadBalancer)
{
  std::string fname = singleton::directories().getLogOutDir() + fileName + ".xml";
  std::ofstream fout;
  if (singleton::mpi().isMainProcessor()) {
    fout.open(fname.c_str(), std::ios::trunc);
    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<XMLContent>\n";
    fout.close();
    fout.clear();
  }

  writeToExistingFile(fname, loadBalancer);

  if (singleton::mpi().isMainProcessor()) {
    fout.open(fname.c_str(), std::ios::app);
    fout << "</XMLContent>\n";
    fout.close();
  }
}


// TODO: Move this method to Cuboid3D<T> class
/// Helper Function to create cuboid parameters for XML tag
template<typename T>
std::string CuboidGeometry3D<T>::_cuboidParameters(Cuboid3D<T> const& cub)
{
  std::stringstream ss;
  ss.flags(std::ios::scientific);
  ss.precision (std::numeric_limits<double>::digits10 + 1);
  ss << " extent=\"";
  for (int i = 0; i<3; i++) {
    ss << cub.getExtend()[i] << " ";
  }

  ss << "\" origin=\"";
  for (int i = 0; i<3; i++) {
    ss << cub.getOrigin()[i] << " ";
  }

  ss << "\" deltaR=\"" << cub.getDeltaR();
  ss << "\" weight=\"" << cub.getWeightValue();
  ss << "\" refinementLevel=\"" << cub.getRefinementLevel() << "\"";
  return ss.str();
}


}  // namespace olb

#endif
