/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2016 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens, Benjamin Förster
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

#ifndef SUPER_LATTICE_INTEGRAL_F_3D_HH
#define SUPER_LATTICE_INTEGRAL_F_3D_HH

#include <cmath>

#include "functors/superLatticeIntegralF3D.h"
#include "functors/indicator/indicatorBaseF3D.hh"
#include "utilities/vectorHelpers.h"
#include "io/ostreamManager.h"

using namespace olb::util;

namespace olb {

template <typename T, typename W>
SuperMin3D<T,W>::SuperMin3D(SuperF3D<T,W>& f,
                            SuperGeometry3D<T>& superGeometry, const int material)
  : SuperF3D<T,W>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Min("+_f.getName()+")";
}

template <typename T, typename W>
bool SuperMin3D<T,W>::operator() (W output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = std::numeric_limits<W>::max();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      for (int iX = 0; iX < nX; ++iX) {
        for (int iY = 0; iY < nY; ++iY) {
          for (int iZ = 0; iZ < nZ; ++iZ) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              W outputTmp[_f.getTargetDim()];
              _f(outputTmp,load.glob(iC),iX,iY,iZ);
              if (outputTmp[i] < output[i]) {
                output[i] = outputTmp[i];
              }
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(output[i], MPI_MIN);
#endif
  }
  return true;
}

template <typename T, typename W>
SuperMax3D<T,W>::SuperMax3D(SuperF3D<T,W>& f,
                            SuperGeometry3D<T>& superGeometry, const int material)
  : SuperF3D<T,W>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Max("+_f.getName()+")";
}

template <typename T, typename W>
bool SuperMax3D<T,W>::operator() (W output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = std::numeric_limits<W>::min();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      for (int iX = 0; iX < nX; ++iX) {
        for (int iY = 0; iY < nY; ++iY) {
          for (int iZ = 0; iZ < nZ; ++iZ) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              W outputTmp[_f.getTargetDim()];
              _f(outputTmp,load.glob(iC),iX,iY,iZ);
              if (outputTmp[i] > output[i]) {
                output[i] = outputTmp[i];
              }
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(output[i], MPI_MAX);
#endif
  }
  return true;
}

template <typename T, typename W>
SuperSum3D<T,W>::SuperSum3D(SuperF3D<T,W>& f,
                            SuperGeometry3D<T>& superGeometry, const int material)
  : SuperF3D<T,W>(f.getSuperStructure(),f.getTargetDim()+1), _f(f),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Sum("+_f.getName()+")";
}

template <typename T, typename W>
bool SuperSum3D<T,W>::operator() (W output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  int numVoxels(0);
  int nX = 0, nY = 0, nZ = 0, iX = 0, iY = 0, iZ = 0;
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = W(0);
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    nX = cGeometry.get(load.glob(iC)).getNx();
    nY = cGeometry.get(load.glob(iC)).getNy();
    nZ = cGeometry.get(load.glob(iC)).getNz();
    for (iX = 0; iX < nX; ++iX) {
      for (iY = 0; iY < nY; ++iY) {
        for (iZ = 0; iZ < nZ; ++iZ) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
            W outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY,iZ);
            for (int i = 0; i < this->getTargetDim()-1; ++i) {
              output[i] += outputTmp[i];
            }
            numVoxels++;
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  output[this->getTargetDim() - 1] = numVoxels;
  return true;
}

template <typename T, typename W>
SuperSumIndicator3D<T,W>::SuperSumIndicator3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, ParticleIndicatorF3D<T,T>& indicator)
  : SuperF3D<T,W>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "Sum("+_f.getName()+")";
}


template <typename T, typename W>
bool SuperSumIndicator3D<T,W>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  int numVoxels(0);
  T outputTmp[_f.getTargetDim()];
  T inside[1];
  Cuboid3D<T>* cub = nullptr;
  int start[3] = {0}, span[3] = {0};

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T(0);
  }

  for (int iC = 0; iC < load.size(); ++iC) {
    int globiC = load.glob(iC);
    cub = &_superGeometry.getCuboidGeometry().get(globiC);

    // check for intersection of cubiod and indicator
    if (cub->getOrigin()[0] <= _indicator.getMax()[0]+_indicator.getPos()[0]
        && cub->getOrigin()[1] <= _indicator.getMax()[1]+_indicator.getPos()[1]
        && cub->getOrigin()[2] <= _indicator.getMax()[2]+_indicator.getPos()[2]
        && _indicator.getMin()[0]+_indicator.getPos()[0] <= cub->getOrigin()[0] + cub->getExtend()[0] * cub->getDeltaR()
        && _indicator.getMin()[1]+_indicator.getPos()[1] <= cub->getOrigin()[1] + cub->getExtend()[1] * cub->getDeltaR()
        && _indicator.getMin()[2]+_indicator.getPos()[2] <= cub->getOrigin()[2] + cub->getExtend()[2] * cub->getDeltaR() ) {

      // compute size of intersection for iteration
      T invDeltaR = 1./cub->getDeltaR();
      for (int k=0; k<3; k++) {
        start[k] = (_indicator.getPos()[k]+_indicator.getMin()[k] - cub->getOrigin()[k]) * invDeltaR;
        if (start[k] < 0) {
          start[k] = 0;
        }
        span[k] = (_indicator.getMax()[k] - _indicator.getMin()[k])*invDeltaR + 3;
        if (span[k] + start[k] > cub->getExtend()[k]) {
          span[k] = cub->getExtend()[k] - start[k];
        }
      }

      for (int iX = start[0]; iX < start[0]+span[0]; iX++) {
        for (int iY = start[1]; iY < start[1]+span[1]; iY++) {
          for (int iZ = start[2]; iZ < start[2]+span[2]; iZ++) {

            _indicator( inside, &(_superGeometry.getPhysR(load.glob(iC), iX, iY, iZ)[0]) );
            if ( !util::nearZero(inside[0]) ) {
              _f(outputTmp,load.glob(iC),iX,iY,iZ);
              for (int i = 0; i < this->getTargetDim()-1; ++i) {
                output[i] += outputTmp[i];
              }
              numVoxels++;
            }
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  output[this->getTargetDim() - 1] = numVoxels;
  return true;
}

template <typename T, typename W>
SuperAverage3D<T,W>::SuperAverage3D(SuperF3D<T,W>& f,
                                    SuperGeometry3D<T>& superGeometry, const int material)
  : SuperF3D<T,double>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Average("+_f.getName()+")";
}


template <typename T, typename W>
bool SuperAverage3D<T,W>::operator() (W output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  int numVoxels(0);
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = W(0);
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();
    for (int iX = 0; iX < nX; ++iX) {
      for (int iY = 0; iY < nY; ++iY) {
        for (int iZ = 0; iZ < nZ; ++iZ) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ)
              == _material) {
            //if (f(load.glob(iC),iX,iY,iZ)[i]!=0) std::cout<< _f(load.glob(iC),iX,iY,iZ)[i] <<std::endl;
            W outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY,iZ);
            for (int i = 0; i < this->getTargetDim()-1 /*f.getTargetDim()*/; ++i) {
              output[i] += outputTmp[i];
            }
            numVoxels++;
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
    output[i] /= numVoxels;
  }
#endif
  output[this->getTargetDim() - 1] = numVoxels;
  return true;
}

template <typename T, typename W>
SuperIntegral3D<T,W>::SuperIntegral3D(SuperF3D<T,W>& f,
                                      SuperGeometry3D<T>& superGeometry, const int material)
  : SuperF3D<T,W>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Integral("+_f.getName()+")";
}

template <typename T, typename W>
bool SuperIntegral3D<T,W>::operator() (W output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = W(0);
  }
  for (int i = 0; i < this->getTargetDim(); ++i) {
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(), 3);
      for (int iX = 0; iX < nX; ++iX) {
        for (int iY = 0; iY < nY; ++iY) {
          for (int iZ = 0; iZ < nZ; ++iZ) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              W outputTmp[_f.getTargetDim()];
              _f(outputTmp,load.glob(iC),iX,iY,iZ);
              output[i] += outputTmp[i]*weight;
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
#endif
  }
  return true;
}

template <typename T, typename W>
SuperLpNorm3D<T,W>::SuperLpNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo,
                                  SuperIndicatorF3D<T>& indicatorF, int p)
  : SuperF3D<T,W>(f.getSuperStructure(),1), _f(f), _superGeometry(geo),
    _indicatorF(indicatorF),  _p(p)
{
  std::ostringstream tmp;
  tmp << _p;
  this->getName() = "L" + tmp.str() + "Norm(" + _f.getName() + ")";
}

template <typename T, typename W>
SuperLpNorm3D<T,W>::SuperLpNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, std::vector<int> materials, int p)
  : SuperLpNorm3D(f, geo, *(new SuperIndicatorMaterial3D<T> (geo, materials)), p)
{};

template <typename T, typename W>
SuperLpNorm3D<T,W>::SuperLpNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, const int material, int p)
  : SuperLpNorm3D(f, geo, std::vector<int> (1, material), p)
{};


template <typename T, typename W>
bool SuperLpNorm3D<T,W>::operator() (W output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  W outputTmp[_f.getTargetDim()];
  bool outputIndicator;

//  int numVoxels(0);
//  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
//    output[iDim] = T(0);
//  }
  output[0] = W(0);
  for (int iC = 0; iC < load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();
    T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(), 3);
    for (int iX = 0; iX < nX; ++iX) {
      for (int iY = 0; iY < nY; ++iY) {
        for (int iZ = 0; iZ < nZ; ++iZ) {
          int pointTmp[4] = { load.glob(iC), iX, iY, iZ };
          _indicatorF(&outputIndicator, pointTmp);
          if (outputIndicator) {
            _f(outputTmp, pointTmp);
            for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
              if     (_p == 0) {
                output[0] = std::max(output[0], std::abs(outputTmp[iDim]));
              } else if (_p == 1) {
                output[0] += std::abs(outputTmp[iDim])*weight;
              } else if (_p == 2) {
                output[0] += outputTmp[iDim]*outputTmp[iDim]*weight;
              } else {
                output[0] += pow(std::abs(outputTmp[iDim]), _p)*weight;
              }
            }
//            numVoxels++;
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  if (_p == 0) {
    singleton::mpi().reduceAndBcast(output[0], MPI_MAX);
  } else {
    singleton::mpi().reduceAndBcast(output[0], MPI_SUM);
  }
//  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
//  output[1] = numVoxels;
  // _p == 1: pass
  // _p == 2: sqrt, else: ^1/_p
  if (_p > 1) {
    output[0] = _p == 2 ? sqrt(output[0]) : pow(output[0], 1. / _p);
  }
  return true;
}
//
//template <typename T, typename W>
//SuperLinfNorm3D<T,W>::SuperLinfNorm3D(SuperF3D<T,W>& f,
//    SuperGeometry3D<T>& superGeometry, const int material)
//  : SuperF3D<T,W>(f.getSuperStructure(),1), _f(f), _superGeometry(superGeometry),
//    _material(material)
//{
//  this->getName() = "LinfNorm("+_f.getName()+")";
//}
//
//template <typename T, typename W>
//bool SuperLinfNorm3D<T,W>::operator() (W output[], const int input[])
//{
//  _f.getSuperStructure().communicate();
//  CuboidGeometry3D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
//  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();
//
//  W outputTmp[_f.getTargetDim()];
//  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
//    output[iDim] = W();
//  }
//
//  for (int iC = 0; iC < load.size(); ++iC) {
//    int nX = cGeometry.get(load.glob(iC)).getNx();
//    int nY = cGeometry.get(load.glob(iC)).getNy();
//    int nZ = cGeometry.get(load.glob(iC)).getNz();
//    for (int iX = 0; iX < nX; ++iX) {
//      for (int iY = 0; iY < nY; ++iY) {
//        for (int iZ = 0; iZ < nZ; ++iZ) {
//          if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ)
//              == _material) {
//            _f(outputTmp, load.glob(iC), iX, iY, iZ);
//            for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
//              if (output[0] < std::abs(outputTmp[iDim]) ) {
//                output[0] = std::abs(outputTmp[iDim]);
//              }
//            }
//          }
//        }
//      }
//    }
//  }
//#ifdef PARALLEL_MODE_MPI
////  int numVoxels(0);
//  singleton::mpi().reduceAndBcast(output[0], MPI_MAX);
////  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
//  //output[1] = numVoxels;
//#endif
//  return true;
//}

template <typename T, typename W>
SuperL1Norm3D<T,W>::SuperL1Norm3D(SuperF3D<T,W>& f,
                                  SuperGeometry3D<T>& geo,
                                  SuperIndicatorF3D<T>& indicatorF)
  : SuperLpNorm3D<T,W>(f, geo, indicatorF, 1)
{
}

template <typename T, typename W>
SuperL1Norm3D<T,W>::SuperL1Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, std::vector<int> materials)
  : SuperLpNorm3D<T,W>(f, geo, materials, 1)
{};

template <typename T, typename W>
SuperL1Norm3D<T,W>::SuperL1Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, const int material)
  : SuperLpNorm3D<T,W>(f, geo, material, 1)
{};


template <typename T, typename W>
SuperL2Norm3D<T,W>::SuperL2Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, SuperIndicatorF3D<T>& indicatorF)
  : SuperLpNorm3D<T,W>(f, geo, indicatorF, 2)
{}

template <typename T, typename W>
SuperL2Norm3D<T,W>::SuperL2Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, std::vector<int> materials)
  : SuperLpNorm3D<T,W>(f, geo, materials, 2)
{};

template <typename T, typename W>
SuperL2Norm3D<T,W>::SuperL2Norm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, const int material)
  : SuperLpNorm3D<T,W>(f, geo, material, 2)
{};


template <typename T, typename W>
SuperLinfNorm3D<T,W>::SuperLinfNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, SuperIndicatorF3D<T>& indicatorF)
  : SuperLpNorm3D<T,W>(f, geo, indicatorF, 0)
{}

template <typename T, typename W>
SuperLinfNorm3D<T,W>::SuperLinfNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, std::vector<int> materials)
  : SuperLpNorm3D<T,W>(f, geo, materials, 0)
{};

template <typename T, typename W>
SuperLinfNorm3D<T,W>::SuperLinfNorm3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& geo, const int material)
  : SuperLpNorm3D<T,W>(f, geo, material, 0)
{};



template<typename T>
SuperGeometryFaces3D<T>::SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry,
    const int material,
    const LBconverter<T>& converter)
  : GenericF<T, int>(7, 0),
    _superGeometry(superGeometry),
    _material(material),
    _converter(converter)
{
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    _blockGeometryFaces.push_back(
      new BlockGeometryFaces3D<T>(_superGeometry.getBlockGeometry(iC),
                                  _material, _converter));
  }
  this->getName() = "superGeometryFaces";
}

template<typename T>
bool SuperGeometryFaces3D<T>::operator()(T output[], const int input[])
{
  _superGeometry.communicate();
  T tmp[7] = { T() };
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim] = T();
  }
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    (*(_blockGeometryFaces[iC]))(tmp, input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += tmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast(output[iDim], MPI_SUM);
  }
#endif
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysDrag3D<T, DESCRIPTOR>::SuperLatticePhysDrag3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry),
    _material(material),
    _faces(_superGeometry, _material, this->_converter),
    _pBoundForce(this->_sLattice, _superGeometry, _material,
                 this->_converter),
    _sumF(_pBoundForce, _superGeometry, _material),
    _factor(
      2.
      / (this->_converter.getCharRho() * this->_converter.getCharU()
         * this->_converter.getCharU()))
{
  this->getName() = "physDrag";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDrag3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  T faces[7] = { 0 };
  T sumF[4] = { 0 };
  _sumF(sumF, input);
  _faces(faces, input);
  //std::cout << faces[0] << std::endl;
  output[0] = _factor * sumF[0] / faces[0];
  output[1] = _factor * sumF[1] / faces[1];
  output[2] = _factor * sumF[2] / faces[2];
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysDragIndicator3D<T, DESCRIPTOR>::SuperLatticePhysDragIndicator3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  ParticleIndicatorSphere3D<T, T>& indicator, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry),
    _indicator(indicator)
{
  this->getName() = "physDrag";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDragIndicator3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  SuperLatticePhysBoundaryForceIndicator3D<T,DESCRIPTOR> pBoundForce(this->_sLattice, _superGeometry, _indicator, this->_converter);
  SuperSumIndicator3D<T,T> sumF(pBoundForce, _superGeometry, _indicator);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());
  T tmp[4] = {};
  sumF(tmp,input);
  output[0] = factor * tmp[0] / _indicator.getDiam();//faces(input)[0];
  output[1] = factor * tmp[1] / _indicator.getDiam();//faces(input)[1];
  output[2] = factor * tmp[2] / _indicator.getDiam();//faces(input)[2];
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysCorrDrag3D<T, DESCRIPTOR>::SuperLatticePhysCorrDrag3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "physCorrDrag";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysCorrDrag3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  SuperGeometryFaces3D<T> faces(_superGeometry, _material, this->_converter);
  SuperLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>  pBoundForce(this->_sLattice, _superGeometry,
      _material, this->_converter);
  SuperSum3D<T,T> sumF(pBoundForce, _superGeometry, _material);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  T sum_tmp[4] = {};
  T face_tmp[7] = {};
  sumF(sum_tmp,input);
  faces(face_tmp,input);
  output[0] = factor * sum_tmp[0] / face_tmp[0];
  output[1] = factor * sum_tmp[1] / face_tmp[1];
  output[2] = factor * sum_tmp[2] / face_tmp[2];
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
  std::vector<T>& u, std::vector<T>& v, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 5),
    _sg(sg),
    _u(u),
    _v(v),
    _origin(A),
    _rad(radius),
    _h(h),
    _analyticalF(f, false)
{
  _mat.push_back(1);
  // normal perpendicular to u and v
  _normal = (T(1) / norm(crossProduct3D(_u, _v))) * crossProduct3D(_u, _v);
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
  std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
  std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 5),
    _sg(sg),
    _u(u),
    _v(v),
    _origin(A),
    _rad(radius),
    _h(h),
    _mat(materials),
    _analyticalF(f, false)
{
  _normal = (T(1) / norm(crossProduct3D(_u, _v))) * crossProduct3D(_u, _v);
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
  std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 5),
    _sg(sg),
    _u(3, 0.),
    _v(3, 0.),
    _origin(A),
    _normal(n),
    _rad(radius),
    _h(h),
    _analyticalF(f, false)
{
  _mat.push_back(1);
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
  std::vector<T>& n, std::vector<T> A, std::list<int> materials, T radius,
  T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 5),
    _sg(sg),
    _u(3, 0.),
    _v(3, 0.),
    _origin(A),
    _normal(n),
    _rad(radius),
    _h(h),
    _mat(materials),
    _analyticalF(f, false)
{
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
  IndicatorCircle3D<T>& circle, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 5),
    _sg(sg),
    _u(3, 0.),
    _v(3, 0.),
    _origin(circle.getCenter()),
    _normal(circle.getNormal()),
    _rad(circle.getRadius()),
    _h(h),
    _vox(0),
    _analyticalF(f, false)
{
  _mat.push_back(1);
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f, SuperGeometry3D<T>& sg,
  IndicatorCircle3D<T>& circle, std::list<int> materials, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 5),
    _sg(sg),
    _u(3, 0.),
    _v(3, 0.),
    _origin(circle.getCenter()),
    _normal(circle.getNormal()),
    _rad(circle.getRadius()),
    _h(h),
    _vox(0),
    _mat(materials),
    _analyticalF(f, false)
{
  init(f);
}

//initialization of member variables
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux3D<T, DESCRIPTOR>::init(
  SuperLatticeF3D<T, DESCRIPTOR>& f)
{
  this->getName() = "SuperLatticeFlux3D";

  //set grid length _h to lattice length
  if (nearZero(_h)) {
    _h = f.getSuperStructure().getCuboidGeometry().getMinDeltaR();
  }

  //define vector _u in plane (perpendicular to _normal)
  if (nearZero(_normal[2])) {
    _u[2] = T(1);
  } else {  //rotation on y-axis by pi/2
    _u[0] = _normal[2];
    _u[2] = -_normal[0];
  }

  //define vector _v in plane perpendicular to _normal and _u
  _v = crossProduct3D(_normal, _u);

  _normal.normalize();  // normalize _normal
  _u.normalize(_h);  // normalize _u and set length _h
  _v.normalize(_h);  // normalize _v and set length _h

  //checking radius
  //maxPhysDist is the diameter of the geometry
  T maxPhysDist = _sg.getStatistics().computeMaxPhysDistance();
  if (_rad < 0 || nearZero(_rad) || _rad > maxPhysDist) {
    _rad = maxPhysDist + _h;
    /*    if (singleton::mpi().getRank() == 0) {
     std::cout << "WARNING: bad radius! Setting radius to maxPhysDist=" << _rad << std::endl;
     }*/
  }
}

//check if point is inside
template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFlux3D<T, DESCRIPTOR>::checkInside(std::vector<T> physR,
    int iC)
{
  std::vector<int> dPos(4, 0);
  //get nearest lattice point
  _sg.getCuboidGeometry().getFloorLatticeR(physR, dPos);
  int iX = dPos[1], iY = dPos[2], iZ = dPos[3];

  //list of material numbers of the eight neighbours of the lattice point
  std::list<int> neighbourCellMaterial;
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY, iZ + 1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY + 1, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY + 1, iZ + 1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX + 1, iY, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX + 1, iY, iZ + 1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX + 1, iY + 1, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX + 1, iY + 1, iZ + 1));

  //if a neighbour has none of the right material numbers it returns false
  bool interpolationPossible = false;
  std::list<int>::iterator i;
  std::list<int>::iterator j;
  for (i = neighbourCellMaterial.begin(); i != neighbourCellMaterial.end();
       ++i) {
    for (j = _mat.begin(); j != _mat.end(); ++j) {
      if (*i == *j) {
        interpolationPossible = true;
      }
    }
    if (interpolationPossible != true) {
      return false;
    } else {
      interpolationPossible = false;
    }
  }
  return true;
}

//summation(integration) of interpolated values
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux3D<T, DESCRIPTOR>::calculate(Vector<T, 3>& flow, int xDir,
    int yDir, int xStart,
    int yStart)
{
  //x and y direction of the plane
  //if positive, start at 0 (default)
  //if negative, start at -1, to prevent double summation at point 0
  if (xDir == -1) {
    xStart = -1;
  }
  if (yDir == -1) {
    yStart = -1;
  }

  //i,j are the affine coordinates of the plane
  int i = xStart;
  int j = yStart;

  //while distance from origin < radius
  while (pow(i * _u[0] + j * _v[0], 2) + pow(i * _u[1] + j * _v[1], 2)
         + pow(i * _u[2] + j * _v[2], 2) < pow(_rad, 2)) {
    while (pow(i * _u[0] + j * _v[0], 2) + pow(i * _u[1] + j * _v[1], 2)
           + pow(i * _u[2] + j * _v[2], 2) < pow(_rad, 2)) {
      T pos[3] = { _origin[0] + i * _u[0] + j * _v[0], _origin[1] + i * _u[1]
                   + j * _v[1], _origin[2] + i * _u[2] + j * _v[2]
                 };
      int iC = _sg.getCuboidGeometry().get_iC(pos[0], pos[1], pos[2], 0);
      if (iC != _sg.getCuboidGeometry().getNc()) {
        if (this->_sg.getLoadBalancer().rank(iC)
            == singleton::mpi().getRank()) {
          //check if point has the right material number
          std::vector<T> vPos(pos, pos + 3);
          if (checkInside(vPos, iC)) {
            T tmp[_analyticalF.getTargetDim()];
            _analyticalF(tmp, pos);
            for (int k = 0; k < 3; k++) {
              if (_analyticalF.getTargetDim() == 3) {  //interpolation
                flow[k] += tmp[k];  //if quantity is three dimensional
              } else {
                flow[k] += tmp[0];  //if quantity is one dimensional
              }
            }
            _vox++;
          }
        }
      }
      j += yDir;
    }
    i += xDir;
    j = yStart;
  }
}

//returns values
template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFlux3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  this->_sLattice.communicate();

  _vox = 0;
  Vector<T, 3> flow;

  //calculate for each of the four quadrants in a plane
  calculate(flow, 1, 1);
  calculate(flow, 1, -1);
  calculate(flow, -1, 1);
  calculate(flow, -1, -1);

  //communicate
#ifdef PARALLEL_MODE_MPI
  for (int j = 0; j < 3; j++) {
    singleton::mpi().reduceAndBcast(flow[j], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(_vox, MPI_SUM);
#endif

  for (int l = 0; l < 3; l++) {
    output[l + 2] = flow[l];      //summation of quantity
    flow[l] *= _h * _h;      //integration
  }

  if (_analyticalF.getTargetDim() == 3) {
    output[0] = flow * _normal;  //projection of flow vector on plane normal (= flux)
  } else {
    output[0] = flow[2];  //summation of quantity, if quantity is one dimensional (= force)
  }

  //number of voxel * grid length^2 (= area)
  output[1] = _vox * _h * _h;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux3D<T, DESCRIPTOR>::print(std::string regionName,
    std::string fluxSiScaleName,
    std::string meanSiScaleName)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v,
  std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _p(sLattice, converter),
    _fluxF(_p, sg, u, v, A, radius, h),
    clout(std::cout, "SuperLatticePhysPressureFlux3D")
{
  this->getName() = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v,
  std::vector<T> A, std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _p(sLattice, converter),
    _fluxF(_p, sg, u, v, A, materials, radius, h),
    clout(std::cout, "SuperLatticePhysPressureFlux3D")
{
  this->getName() = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _p(sLattice, converter),
    _fluxF(_p, sg, n, A, radius, h),
    clout(std::cout, "SuperLatticePhysPressureFlux3D")
{
  this->getName() = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A,
  std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _p(sLattice, converter),
    _fluxF(_p, sg, n, A, materials, radius, h),
    clout(std::cout, "SuperLatticePhysPressureFlux3D")
{
  this->getName() = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<T>& circle, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _p(sLattice, converter),
    _fluxF(_p, sg, circle, h),
    clout(std::cout, "SuperLatticePhysPressureFlux3D")
{
  this->getName() = "SuperLatticePhysPressureFlux3D";
}
template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<T>& circle,
  std::list<int> materials, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _p(sLattice, converter),
    _fluxF(_p, sg, circle, materials, h),
    clout(std::cout, "SuperLatticePhysPressureFlux3D")
{
  this->getName() = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  _fluxF(output, input);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::print(
  std::string regionName, std::string fluxSiScaleName,
  std::string meanSiScaleName)
{
  int input[1] = { };
  T output[_fluxF.getTargetDim()];
  this->operator()(output, input);
  if (regionName != "") {
    clout << "regionName=" << regionName << "; regionSize[m^2]=" << output[1]
          << std::flush;
  } else {
    clout << "regionSize[m^2]=" << output[1] << std::flush;
  }
  if (singleton::mpi().isMainProcessor()) {
    if (fluxSiScaleName == "MN") {
      std::cout << "; force[MN]=" << output[0] / T(1.e6) << std::flush;
    } else if (fluxSiScaleName == "kN") {
      std::cout << "; force[kN]=" << output[0] / T(1.e3) << std::flush;
    } else {
      std::cout << "; force[N]=" << output[0] << std::flush;
    }
    if ( meanSiScaleName == "mmHg" ) {
      std::cout << "; meanPressure[mmHg]=" << std::abs(output[0])/output[1]/T(133.322) << std::endl;
    } else {
      std::cout << "; meanPressure[Pa]=" << std::abs(output[0])/output[1] << std::endl;
    }
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v,
  std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _vel(sLattice, converter),
    _fluxF(_vel, sg, u, v, A, radius, h),
    clout(std::cout, "SuperLatticePhysVelocityFlux3D")
{
  this->getName() = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v,
  std::vector<T> A, std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _vel(sLattice, converter),
    _fluxF(_vel, sg, u, v, A, materials, radius, h),
    clout(std::cout, "SuperLatticePhysVelocityFlux3D")
{
  this->getName() = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, radius, h),
    clout(std::cout, "SuperLatticePhysVelocityFlux3D")
{
  this->getName() = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A,
  std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, materials, radius, h),
    clout(std::cout, "SuperLatticePhysVelocityFlux3D")
{
  this->getName() = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<T>& circle, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _vel(sLattice, converter),
    _fluxF(_vel, sg, circle, h),
    clout(std::cout, "SuperLatticePhysVelocityFlux3D")
{
  this->getName() = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<T>& circle,
  std::list<int> materials, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5),
    _vel(sLattice, converter),
    _fluxF(_vel, sg, circle, materials, h),
    clout(std::cout, "SuperLatticePhysVelocityFlux3D")
{
  this->getName() = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  _fluxF(output, input);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::print(
  std::string regionName, std::string fluxSiScaleName,
  std::string meanSiScaleName)
{
  int input[1] = { };
  T output[_fluxF.getTargetDim()];
  this->operator()(output, input);
  if (regionName != "") {
    clout << "regionName=" << regionName << "; regionSize[m^2]=" << output[1]
          << std::flush;
  } else {
    clout << "regionSize[m^2]=" << output[1] << std::flush;
  }
  if (singleton::mpi().isMainProcessor()) {
    if (fluxSiScaleName == "ml/s") {
      std::cout << "; volumetricFlowRate[ml/s]=" << output[0] * T(1.e6)
                << std::flush;
    } else if (fluxSiScaleName == "l/s") {
      std::cout << "; volumetricFlowRate[l/s]=" << output[0] * T(1.e3)
                << std::flush;
    } else {
      std::cout << "; volumetricFlowRate[m^3/s]=" << output[0] << std::flush;
    }
    if (meanSiScaleName == "mm/s") {
      std::cout << "; meanVelocity[mm/s]=" << output[0] / output[1] * T(1.e3)
                << std::endl;
    } else {
      std::cout << "; meanVelocity[m/s]=" << output[0] / output[1] << std::endl;
    }
  }
}

}  // end namespace olb

#endif
