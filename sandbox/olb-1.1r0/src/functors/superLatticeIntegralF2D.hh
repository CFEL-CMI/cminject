/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause,
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

#ifndef SUPER_LATTICE_INTEGRAL_F_2D_HH
#define SUPER_LATTICE_INTEGRAL_F_2D_HH

#include<vector>
#include<cmath>

#include "functors/superLatticeIntegralF2D.h"
#include "functors/blockLatticeIntegralF2D.h"
#include "utilities/vectorHelpers.h"
#include "core/vector.h"

using namespace olb::util;

namespace olb {


template <typename T>
SuperMax2D<T>::SuperMax2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Max("+_f.getName()+")";
}

template <typename T>
bool SuperMax2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            T outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY);
            if (fabs(outputTmp[i]) > output[i]) {
              output[i] = fabs(outputTmp[i]);
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

template <typename T>
SuperMin2D<T>::SuperMin2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Min("+_f.getName()+")";
}

template <typename T>
bool SuperMin2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            T outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY);
            if (fabs(outputTmp[i]) < output[i]) {
              output[i] = fabs(outputTmp[i]);
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


template <typename T>
SuperSum2D<T>::SuperSum2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Sum("+_f.getName()+")";
}

template <typename T>
bool SuperSum2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  int numVoxels(0);
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
  }
  for (int iC=0; iC<load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          T outputTmp[_f.getTargetDim()];
          _f(outputTmp,load.glob(iC),iX,iY);
          for (int i = 0; i < this->getTargetDim()-1 /*f.getTargetDim()*/; ++i) {
            output[i] += outputTmp[i];
          }
          numVoxels++;
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
  output[this->getTargetDim()-1] = numVoxels;
  return true;
}


template <typename T>
SuperSumIndicator2D<T>::SuperSumIndicator2D(SuperF2D<T>& f,
    SuperGeometry2D<T>& superGeometry, ParticleIndicatorF2D<T,T>& indicator)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "Sum("+_f.getName()+")";
}

template <typename T>
bool SuperSumIndicator2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = 0.;
  }

  T physR[2];
  T inside[1];
  int numVoxels(0);
  T outputTmp[_f.getTargetDim()];
  Cuboid2D<T>* cub = nullptr;
  int start[2] = {0}, span[2] = {0};
  for (int iC = 0; iC < load.size(); ++iC) {
    int globiC = load.glob(iC);
    cub = &_superGeometry.getCuboidGeometry().get(globiC);
    if (! (cub->get_globPosX() > _indicator.getPos()[0]+_indicator.getMax()[0] ||
           cub->get_globPosY() > _indicator.getPos()[1]+_indicator.getMax()[1] ||
           _indicator.getPos()[0]+_indicator.getMin()[0] > cub->get_globPosX() + cub->getExtend()[0] * cub->getDeltaR() ||
           _indicator.getPos()[1]+_indicator.getMin()[1] > cub->get_globPosY() + cub->getExtend()[1] * cub->getDeltaR())) {
      for (int k=0; k<2; k++) {
        start[k] = (_indicator.getPos()[k]+_indicator.getMin()[k] - cub->getOrigin()[k]) /
                   cub->getDeltaR();
        if (start[k] < 0) {
          start[k] = 0;
        }
        span[k] = (_indicator.getMax()[k] -
                   _indicator.getMin()[k])/cub->getDeltaR() + 3;
        if (span[k] + start[k] > cub->getExtend()[k]) {
          span[k] = cub->getExtend()[k] - start[k];
        }
      }
      for (int iX = start[0]; iX < start[0]+span[0]; iX++) {
        for (int iY = start[1]; iY < start[1]+span[1]; iY++) {
          if (_superGeometry.get(globiC, iX, iY) == 1) {
            cub->getPhysR(physR,iX,iY);
            _indicator(inside, physR);
            if ( !util::nearZero(inside[0]) ) {
              _f(outputTmp,globiC,iX,iY);
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
  output[this->getTargetDim()-1] = numVoxels;
  return true;
}


template <typename T>
SuperIntegral2D<T>::SuperIntegral2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()), _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Integral("+_f.getName()+")";
}

template <typename T>
bool SuperIntegral2D<T>::operator() (T output[], const int input[])
{
  //  f.getSuperStructure().communicate();
  //  CuboidGeometry2D<T>& cGeometry = f.getSuperStructure().getCuboidGeometry();
  //  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();

  //  std::vector<T> tmp(this->_n, T() );
  //  for (int i=0; i<this->_n; ++i) {
  //    for (int iC=0; iC<load.size(); ++iC) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  ////      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      T weight = pow(this->superGeometry.getDeltaR(),3);
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  ////          for (int iZ=0; iZ<nZ; iZ++) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  ////            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  ////            if (this->superGeometry.getMaterial(globX, globY) == material) {
  ////              tmp[i]+=f(load.glob(iC),iX,iY)[i]*weight;
  ////            }
  ////            if (this->superGeometry.getMaterial(globX, globY, globZ) == material) {
  ////              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*weight;
  ////            }

  ////          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  //  return tmp;
  return false;
}

template <typename T>
SuperL1Norm2D<T>::SuperL1Norm2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),1), _f(f), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "L1("+_f.getName()+")";
}

template <typename T>
bool SuperL1Norm2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  T outputTmp[_f.getTargetDim()];
//  int numVoxels(0);
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = T();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),2);
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          _f(outputTmp,load.glob(iC),iX,iY);
          for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
            output[0] += fabs(outputTmp[iDim])*weight;
          }
//          numVoxels++;
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(output[0], MPI_SUM);
//  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  //output[1] = numVoxels;
  return true;
}


template <typename T>
SuperL2Norm2D<T>::SuperL2Norm2D(SuperF2D<T>& f,
                                SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),1), _f(f), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "L2Norm("+_f.getName()+")";
}

template <typename T>
bool SuperL2Norm2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  T outputTmp[_f.getTargetDim()];
//  int numVoxels(0);
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = T();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),2);
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          _f(outputTmp,load.glob(iC),iX,iY);
          for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
            output[0] += outputTmp[iDim]*outputTmp[iDim]*weight;
          }
//          numVoxels++;
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(output[0], MPI_SUM);
//  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  //output[1] = numVoxels;
  output[0] = sqrt(output[0]);
  return true;
}


template <typename T>
SuperLinfNorm2D<T>::SuperLinfNorm2D(SuperF2D<T>& f,
                                    SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),1), _f(f), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "LinfNorm("+_f.getName()+")";
}

template <typename T>
bool SuperLinfNorm2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  T outputTmp[_f.getTargetDim()];
//  int numVoxels(0);
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = T();
  }

  for (int iC = 0; iC < load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          _f(outputTmp,load.glob(iC),iX,iY);
          for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
            if (fabs(output[0]) < fabs(outputTmp[iDim]) ) {
              output[0] = fabs(outputTmp[iDim]);
            }
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(output[0], MPI_MAX);
//  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  //output[1] = numVoxels;
  return true;
}

template <typename T>
SuperL222D<T>::SuperL222D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "L22("+_f.getName()+")";
}

template <typename T>
bool SuperL222D<T>::operator() (T output[], const int input[])
{
  //  f.getSuperStructure().communicate();
  //  CuboidGeometry2D<T>& cGeometry = f.getSuperStructure().getCuboidGeometry();
  //  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();

  //  std::vector<T> tmp(this->_n, T() );
  //  for (int i=0; i<this->_n; ++i) {

  //    for (int iC=0; iC<load.size(); ++iC) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  ////      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      T weight = pow(this->superGeometry.getDeltaR(),3);
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  ////          for (int iZ=0; iZ<nZ; iZ++) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  ////            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  ////            if (this->superGeometry.getMaterial(globX, globY) == material) {
  ////              //std::cout << f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight << std::endl;
  ////              tmp[i]+=f(load.glob(iC),iX,iY)[i]*f(load.glob(iC),iX,iY)[i]*weight;
  ////            }
  ////            if (this->superGeometry.getMaterial(globX, globY, globZ) == material) {
  ////              //std::cout << f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight << std::endl;
  ////              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight;
  ////            }

  ////          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  //  return tmp;
  return false;
}



template <typename T>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry,
    const int material, const LBconverter<T>& converter)
  : GenericF<T,int>(7,3), _superGeometry(superGeometry), _material(material),
    _converter(converter)
{
  this->getName() = "superGeometryFaces";
}

template <typename T>
bool SuperGeometryFaces2D<T>::operator() (T output[], const int input[])
{
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim]=T();
  }
  _superGeometry.communicate();
  std::vector<T> counter(7,T());
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    BlockGeometryFaces2D<T> f(_superGeometry.getBlockGeometry(iC), _material, _converter);
    T outputTmp[f.getTargetDim()];
    f(outputTmp,input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += outputTmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast( output[iDim], MPI_SUM);
  }
#endif
  return true;
}


template <typename T>
SuperGeometryFacesIndicator2D<T>::SuperGeometryFacesIndicator2D(SuperGeometry2D<T>& superGeometry,
    SmoothIndicatorCircle2D<T,T>& indicator, const int material, const LBconverter<T>& converter)
  : GenericF<T,int>(7,3), _superGeometry(superGeometry), _indicator(indicator), _material(material),
    _converter(converter)
{
  this->getName() = "superGeometryFacesInd";
}

template <typename T>
bool SuperGeometryFacesIndicator2D<T>::operator() (T output[], const int input[])
{
  _superGeometry.communicate();
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim]=T();
  }
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    BlockGeometryFacesIndicator2D<T> f(_superGeometry.getBlockGeometry(iC), _indicator, _material, _converter);
    T outputTmp[f.getTargetDim()];
    f(outputTmp,input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += outputTmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast( output[iDim], MPI_SUM);
  }
#endif
  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDrag2D<T,DESCRIPTOR>::SuperLatticePhysDrag2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physDrag";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  SuperGeometryFaces2D<T> faces(_superGeometry, _material, this->_converter);
  SuperLatticePhysBoundaryForce2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _material, this->_converter);
  SuperSum2D<T> sumF(f, _superGeometry, _material);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  T facesTmp[faces.getTargetDim()], sumFTmp[sumF.getTargetDim()];
  sumF(sumFTmp, input);
  faces(facesTmp, input);
  output[0] = factor * sumFTmp[0] / facesTmp[0];
  output[1] = factor * sumFTmp[1] / facesTmp[1];

  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDragIndicator2D<T,DESCRIPTOR>::SuperLatticePhysDragIndicator2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 ParticleIndicatorF2D<T,T>& indicator, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physDragIndicator";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysDragIndicator2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  SuperLatticePhysBoundaryForceIndicator2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _indicator, this->_converter);
  SuperSumIndicator2D<T> sumF(f, _superGeometry, _indicator);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  T sumFTmp[sumF.getTargetDim()];
  sumF(sumFTmp, input);
  output[0] = factor * sumFTmp[0] / _indicator.getDiam();
  output[1] = factor * sumFTmp[1] / _indicator.getDiam();

  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDragIndicator2D_2<T,DESCRIPTOR>::SuperLatticePhysDragIndicator2D_2
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 SmoothIndicatorF2D<T,T>& indicator, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physDragIndicator";
}

/*
template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysDragIndicator2D_2<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  SuperLatticePhysBoundaryForceIndicator2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _indicator, this->_converter);
//  SuperLatticePhysVolumeForceIndicator2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _indicator, this->_converter);
  SuperSumIndicator2D<T> sumF(f, _superGeometry, _indicator);

  T sumFTmp[sumF.getTargetDim()];
  sumF(output, input);

  return true;
}
*/

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::SuperLatticePhysCorrDrag2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physCorrDrag";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  SuperGeometryFaces2D<T> faces(superGeometry, material, this->converter);

  //  SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR> f(this->sLattice, superGeometry, material, this->converter);
  //  SuperSum2D<T> sumF(f, superGeometry, material);

  //  T factor = 2. / (this->converter.getCharRho() * this->converter.getCharU() * this->converter.getCharU());

  //  std::vector<T> drag(2,T());
  //  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  //  drag[1] = factor * sumF(input)[1] / faces(input)[1];
  ////  drag[2] = factor * sumF(input)[2] / faces(input)[2];
  //  return drag;
  return false;
}



template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux2D<T, DESCRIPTOR>::SuperLatticeFlux2D(SuperLatticeF2D<T, DESCRIPTOR>& f,
    SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(f.getSuperLattice(),4), _sg(sg), _u(2, 0.),
    _origin(A), _normal(n), _rad(radius), _h(h), _vox(0), _analyticalF(f)
{
  init(f);
  _mat.push_back(1);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux2D<T, DESCRIPTOR>::SuperLatticeFlux2D(SuperLatticeF2D<T, DESCRIPTOR>& f,
    SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, std::list<int> materials,
    T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(f.getSuperLattice(),4), _sg(sg), _u(2, 0.),
    _origin(A), _normal(n), _rad(radius), _h(h), _vox(0), _mat(materials), _analyticalF(f)
{
  init(f);
}


//initialization of member variables
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux2D<T, DESCRIPTOR>::init(SuperLatticeF2D<T, DESCRIPTOR>& f)
{
  if (nearZero(_h) ) {
    _h = f.getSuperLattice().getCuboidGeometry().getMinDeltaR();
  }
  //normalize normal
  _normal.normalize();

  //rotation of normal 90° clockwise
  _u[0] = _normal[1];
  _u[1] = -_normal[0];
  T tmp = _h/norm(Vector<T,2>(_u));
  _u[0] *= tmp;
  _u[1] *= tmp;

  //checking radius
  //maxPhysDist is the diameter of the geometry
  T maxPhysDist = pow(_sg.getStatistics().getMaxPhysR(1)[0]-_sg.getStatistics().getMinPhysR(1)[0],2)
                  +pow(_sg.getStatistics().getMaxPhysR(1)[1]-_sg.getStatistics().getMinPhysR(1)[1],2);
  maxPhysDist = sqrt(maxPhysDist);
  if (_rad < 0 || nearZero(_rad) || _rad > maxPhysDist) {
    _rad = maxPhysDist + _h;
    if (singleton::mpi().getRank() == 0) {
      std::cout << "WARNING: bad radius! Setting radius to maxPhysDist=" << _rad << std::endl;
    }
  }
}


//check if point is inside
template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFlux2D<T, DESCRIPTOR>::checkInside(std::vector<T> physR, int iC)
{
  std::vector<int> dPos(4,0);
  _sg.getCuboidGeometry().getFloorLatticeR(physR, dPos);
  int iX = dPos[1], iY = dPos[2];

  //list of material numbers of the four neighbours
  std::list<int> neighbourCellMaterial;
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY+1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY+1));

  //if a neighbour has none of the right material numbers of _mat it returns false
  bool interpolationPossible = false;
  std::list<int>::iterator i;
  std::list<int>::iterator j;
  for (i = neighbourCellMaterial.begin(); i != neighbourCellMaterial.end(); ++i) {
    for (j = _mat.begin(); j != _mat.end(); ++j) {
      if (*i == *j) {
        interpolationPossible = true;
      }
    }
    if (interpolationPossible!=true) {
      return false;
    } else {
      interpolationPossible=false;
    }
  }
  return true;
}

//summation(integration) of interpolated values
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux2D<T, DESCRIPTOR>::calculate(std::vector<T>& flow, int xDir,
    int xStart)
{
  //xDir is direction of the line
  //if positive, start at 0 (default)
  //if negative, start at -1, to prevent double summation at point 0
  if (xDir == -1) {
    xStart = -1;
  }

  int i = xStart;

  //while distance from origin < radius
  while (pow(i*_u[0],2)+pow(i*_u[1],2) < pow(_rad,2)) {
    T pos[2] = {_origin[0] + i*_u[0], _origin[1] + i*_u[1]};
    int iC = _sg.getCuboidGeometry().get_iC(pos[0],pos[1], 0);
    if (this->_sg.getLoadBalancer().rank(iC) == singleton::mpi().getRank()) {
      //check if point has the right material number
      std::vector<T> vPos(pos,pos+2);
      if (checkInside(vPos,iC)) {
        T tmp[_analyticalF.getTargetDim()];
        _analyticalF(tmp,pos);
        for (int k = 0; k < 2; k++) {
          if (_analyticalF.getTargetDim() == 2) {
            flow[k] += tmp[k];  //interpolation
          } else {
            flow[k] += tmp[0];
          }
        }
        _vox++;
      }
    }
    i += xDir;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFlux2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  this->_sLattice.communicate();

  _vox = 0;
  std::vector<T> flow(4, T());
  std::vector<T> tmpFlow(4, T());
  int iXp;
  int iXn;
  iXp = 1;
  iXn = -1;

  //calculate for both direction of the line
  calculate(flow, iXp);
  calculate(flow, iXn);

  //communicate
#ifdef PARALLEL_MODE_MPI
  for (int j = 0; j < 2; j++) {
    singleton::mpi().reduceAndBcast(flow[j], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(_vox, MPI_SUM);
#endif

  //quadrature
  for (int l = 0; l < 2; l++) {
    output[l+2] = flow[l]; //summation of quantity
    flow[l] *= _h;          //integration
  }

  if (_analyticalF.getTargetDim() == 2) {
    output[0] = flow[0]*_normal[0] + flow[1]*_normal[1];  //projection of flow vector on plane normal (= flux)
  } else {
    output[0] = flow[1];  //summation of quantity, if quantity is a scalar (= force)
  }

  //number of voxel * grid length (= length)
  output[1] = _vox*_h;

  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux2D<T, DESCRIPTOR>::print(std::string regionName,
    std::string fluxSiScaleName, std::string meanSiScaleName ) { }



template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux2D
(SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
 SuperGeometry2D<T>& sg, std::vector<T>& u, std::vector<T> A, T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _p(sLattice, converter),
    _fluxF(_p, sg, u, A, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux2D")
{
  this->getName() = "SuperLatticePhysPressureFlux2D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux2D
(SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
 SuperGeometry2D<T>& sg, std::vector<T>& u, std::vector<T> A, std::list<int> materials,
 T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _p(sLattice, converter),
    _fluxF(_p, sg, u, A, materials, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux2D")
{
  this->getName() = "SuperLatticePhysPressureFlux2D";
}


template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  _fluxF(output,input);
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::print(std::string regionName,
    std::string fluxSiScaleName, std::string meanSiScaleName)
{
  int input[1]= {};
  T output[_fluxF.getTargetDim()];
  this->operator()(output,input);
  if ( regionName != "") {
    clout << "regionName=" << regionName << "; regionSize[m]=" << output[1] << std::flush;
  } else {
    clout << "regionSize[m]=" << output[1] << std::flush;
  }
  if (singleton::mpi().isMainProcessor() ) {
    if ( fluxSiScaleName == "MN" ) {
      std::cout << "; force[MN]=" << output[0]/T(1.e6) << std::flush;
    } else if ( fluxSiScaleName == "kN") {
      std::cout << "; force[kN]=" << output[0]/T(1.e3) << std::flush;
    } else {
      std::cout << "; force[N]=" << output[0] << std::flush;
    }
    if ( meanSiScaleName == "mmHg" ) {
      std::cout << "; meanPressure[mmHg]=" << fabs(output[0])/output[1]/T(133.322) << std::endl;
    } else {
      std::cout << "; meanPressure[Pa]=" << fabs(output[0])/output[1] << std::endl;
    }
  }
}



template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux2D
(SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
 SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux2D")
{
  this->getName() = "SuperLatticePhysVelocityFlux2D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux2D
(SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
 SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, std::list<int> materials,
 T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, materials, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux2D")
{
  this->getName() = "SuperLatticePhysVelocityFlux2D";
}


template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  _fluxF(output,input);
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::print(std::string regionName,
    std::string fluxSiScaleName, std::string meanSiScaleName )
{
  int input[1] = {};
  T output[_fluxF.getTargetDim()];
  this->operator()(output,input);
  if ( regionName != "") {
    clout << "regionName=" << regionName << "; regionSize[m]=" << output[1] << std::flush;
  } else {
    clout << "regionSize[m]=" << output[1] << std::flush;
  }
  if (singleton::mpi().isMainProcessor() ) {
    std::cout << "; flowRate[m^2/s]=" << output[0] << std::flush;
    if ( meanSiScaleName == "mm/s" ) {
      std::cout << "; meanVelocity[mm/s]=" << output[0]/output[1]*T(1.e3) << std::endl;
    } else {
      std::cout << "; meanVelocity[m/s]=" << output[0]/output[1] << std::endl;
    }
  }
}



#endif

} //end namespace olb
