/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens, Benjamin Förster, Marie-Luise Maier
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

#ifndef INTERPOLATION_F_3D_HH
#define INTERPOLATION_F_3D_HH

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "functors/interpolationF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include <algorithm>

namespace olb {

/// a class used to convert a block functor to an analytical functor
template <typename T, typename W>
AnalyticalFfromBlockF3D<T,W>::AnalyticalFfromBlockF3D(BlockF3D<W>& f,
    Cuboid3D<T>& cuboid)
  : AnalyticalF3D<T,W>(f.getTargetDim()), _f(f), _cuboid(cuboid)
{
  this->getName() = "fromBlockF";
}


template <typename T, typename W>
bool AnalyticalFfromBlockF3D<T,W>::operator()(W output[], const T physC[])
{
  int latticeR[3];
  _cuboid.getFloorLatticeR(latticeR, physC);

  Vector<T,3> physRiC;
  Vector<T,3> physCv (physC);
  Vector<W,3> d, e;
  W output_tmp[3];

  _cuboid.getPhysR(physRiC.data, latticeR[0], latticeR[1], latticeR[2]);
  T dr = 1/_cuboid.getDeltaR();

  // compute weights
  d = (physCv - physRiC) * dr;
  e = 1. - d;

  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = W();
    output_tmp[iD] = W();
  }

  //0=1=2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*e[2];
  }

  latticeR[1]++;
  //0=1+2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*e[2];
  }

  latticeR[0]++;
  latticeR[1]--;
  //0+1=2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*e[2];
  }

  latticeR[1]++;
  //0+1+2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*e[2];
  }

  latticeR[0]--;
  latticeR[1]--;
  latticeR[2]++;
  //0=1=2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*d[2];
  }

  latticeR[1]++;
  //0=1+2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*d[2];
  }

  latticeR[0]++;
  latticeR[1]--;
  //0+1=2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*d[2];
  }

  latticeR[1]++;
  //0+1+2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*d[2];
  }

  return true;
}



template <typename T, typename W>
SpecialAnalyticalFfromBlockF3D<T,W>::SpecialAnalyticalFfromBlockF3D(
  BlockF3D<W>& f, Cuboid3D<T>& cuboid,
  Vector<T,3> delta)
  : AnalyticalF3D<T,W>(f.getTargetDim()), _f(f), _cuboid(cuboid), _delta(delta)
{
  this->getName() = "fromBlockF";
}


template <typename T, typename W>
bool SpecialAnalyticalFfromBlockF3D<T,W>::operator()(W output[],
    const T physC[])
{
  Vector<T,3> origin = _cuboid.getOrigin();

  // scale physC in all 3 dimensions
  Vector<T,3> physCv;
  for (int i=0; i<3; i++) {
    physCv[i] = origin[i] + (physC[i] - origin[i]) * ( _cuboid.getDeltaR() /
                _delta[i] );
  }


  int latticeR[3];
  // NOTE: latticeR[i] may be < 0 for some i! -> calculate manually below
  //  _cuboid.getFloorLatticeR(latticeR, physC);
  // std::cout << singleton::mpi().getRank() << "] physC[2] = " << physC[2]
  // << " - origin[2] = " << origin[2] << " - LatticeR = " << latticeR[0] <<
  // " : " << latticeR[1] << " : " << latticeR[2] << std::endl;
//  Vector<T,3> translatedPhysC;
  for (int i=0; i<3; i++) {
//    if( _cuboid.getDeltaR() / _delta[i]  == T(0) ) {
//      std::cout << "00000000000000000000" << std::endl;
//    }
//    translatedPhysC[i] = _cuboid.getDeltaR() / _delta[i] * physC[i];
    latticeR[i] = std::max((int)floor( (physCv[i] - origin[i])/
                                       _cuboid.getDeltaR()), 0);
//    latticeR[i] = std::max((int)floor( (physC[i] - origin[i]) * _delta[i] /
    //(_cuboid.getDeltaR() * _cuboid.getDeltaR())), 0);
//    latticeR[i] = (int)floor( (physC[i] - origin[i]) * _delta[i] /
    //(_cuboid.getDeltaR() * _cuboid.getDeltaR()));
  }
//  if(latticeR[2]) {
////    std::cout << "origin[2] = " << origin[2] << std::endl;
////    std::cout << "_cuboid.getDeltaR() / _delta[2] = " <<
  // _cuboid.getDeltaR() / _delta[2] << std::endl;
////    std::cout << "_cuboid.getDeltaR() / _delta[2] * physC[2] = " <<
  // _cuboid.getDeltaR() / _delta[2] * physC[2] << std::endl;
////    std::cout << "(translatedPhysC[2] - origin[2])/_cuboid.getDeltaR() =
  // " << (translatedPhysC[2] - origin[2])/_cuboid.getDeltaR() << std::endl;
//    std::cout << "latticeR[" << 2 << "] = " << latticeR[2] << std::endl;
//  }

  Vector<T,3> physRiC;
  Vector<W,3> d, e;
  W output_tmp[3];
//  _f(output, latticeR);
//  return true;
  Vector<T,3> latticeRv;
  for (int i=0; i<3; i++) {
//    dr[i] = 1./_delta[i];
    latticeRv[i] = (T) latticeR[i];
  }

//  physRiC = origin + latticeRv * _delta;
  physRiC = origin + latticeRv * _cuboid.getDeltaR();

//  _cuboid.getPhysR(physRiC.data, latticeR[0], latticeR[1], latticeR[2]);



//  Vector<T,3> dr;
//  for(int i=0;i<3;i++){
////    dr[i] = 1./_delta[i];
//    dr[i] = 1./_cuboid.getDeltaR();
//  }
  T dr = 1. / _cuboid.getDeltaR();

  // compute weights
  d = (physCv - physRiC) * dr;
  e = 1. - d;

  // compute weights
  d = (physCv - physRiC) * dr;
  e = 1. - d;

  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = W();
    output_tmp[iD] = W();
  }

  //0=1=2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*e[2];
  }

  latticeR[1]++;
  //0=1+2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*e[2];
  }

  latticeR[0]++;
  latticeR[1]--;
  //0+1=2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*e[2];
  }

  latticeR[1]++;
  //0+1+2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*e[2];
  }

  latticeR[0]--;
  latticeR[1]--;
  latticeR[2]++;
  //0=1=2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*d[2];
  }

  latticeR[1]++;
  //0=1+2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*d[2];
  }

  latticeR[0]++;
  latticeR[1]--;
  //0+1=2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*d[2];
  }

  latticeR[1]++;
  //0+1+2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*d[2];
  }

  return true;
}


/// a class used to convert lattice functions to analytical functions
template <typename T, typename W>
AnalyticalFfromSuperF3D<T,W>::AnalyticalFfromSuperF3D(SuperF3D<T,W>& f,
    bool communicateToAll, int overlap,
    bool communicateOverlap)
  : AnalyticalF3D<T,W>(f.getTargetDim()), _f(f),
    _cuboidGeometry(f.getSuperStructure().getCuboidGeometry()),
    _communicateToAll(communicateToAll), _overlap(overlap),
    _communicateOverlap(communicateOverlap)
{
  this->getName() = "fromSuperF";
  //std::cout << _f.getName() << std::endl;
  if (overlap == -1) {
    _overlap = f.getSuperStructure().getOverlap();
  }
  /*if (&(_f.getBlockF(0))!=NULL ) {
     for (int iC = 0; iC < _f.getSuperStructure().getLoadBalancer().size();
     iC++ ) {
       int iCglob = _f.getSuperStructure().getLoadBalancer().glob(iC);
       this->_analyticalFfromBlockF.push_back(new
       AnalyticalFfromBlockF3D<T>(_f.getBlockF(iC), _f.getSuperStructure().
       getCuboidGeometry().get(iCglob) ) );
     }
   }*/
}


template <typename T, typename W>
AnalyticalFfromSuperF3D<T,W>::~AnalyticalFfromSuperF3D()
{
  if (_analyticalFfromBlockF.size() != 0) {
    for (unsigned int iC = 0; iC < _analyticalFfromBlockF.size(); ++iC) {
      delete _analyticalFfromBlockF[iC];
    }
    _analyticalFfromBlockF.resize(0);
  }
}

template <typename T, typename W>
bool AnalyticalFfromSuperF3D<T,W>::operator()(W output[], const T physC[])
{
  int latticeR[4];
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = W();
  }
  if (!(_cuboidGeometry.getLatticeR(latticeR, physC))) {
    return false;
  }
  // convert to lattice coordinates
  W d[3];

  if (_communicateOverlap) {
    _f.getSuperStructure().communicate();
  }

  int locX, locY, locZ;

  int dataSize = 0;
  int dataFound = 0;

  int latticeC[4] = {};

  for (int iC = 0; iC < _f.getSuperStructure().getLoadBalancer().size(); ++iC) {
    latticeC[0] = _f.getSuperStructure().getLoadBalancer().glob(iC);
    _cuboidGeometry.get(latticeC[0]).getFloorLatticeR(latticeR, physC);
    if (latticeR[0] >= -_overlap && latticeR[0] + 1 <
        _cuboidGeometry.get(latticeC[0]).getNx() + _overlap &&
        latticeR[1] >= -_overlap && latticeR[1] + 1 <
        _cuboidGeometry.get(latticeC[0]).getNy() + _overlap &&
        latticeR[2] >= -_overlap && latticeR[2] + 1 <
        _cuboidGeometry.get(latticeC[0]).getNz() + _overlap ) {
      if (_analyticalFfromBlockF.size() != 0 ) {
        (*(_analyticalFfromBlockF[iC]))(output, physC);
      } else {
        locX = latticeR[0];
        locY = latticeR[1];
        locZ = latticeR[2];

        T physRiC[3];
        _cuboidGeometry.get(latticeC[0]).getPhysR(physRiC, locX, locY, locZ);

        d[0] = (physC[0] - physRiC[0]);
        d[1] = (physC[1] - physRiC[1]);
        d[2] = (physC[2] - physRiC[2]);

        d[0] /= _cuboidGeometry.get(latticeC[0]).getDeltaR();
        d[1] /= _cuboidGeometry.get(latticeC[0]).getDeltaR();
        d[2] /= _cuboidGeometry.get(latticeC[0]).getDeltaR();

        //if (d[0]*d[0]>1||d[1]*d[1]>1||d[2]*d[2]>1)
        //cout << d[0] << " " << d[1] << " " << d[2] << endl;
        //if (locX<0||locY<0||locZ<0)
        //cout << locX << " " << locY << " " << locZ << "_overlap="
        //<<_overlap<< endl;

        //      std::vector<T> output(f(latticeC).size(), T());
        //std::vector<T> output_tmp(_f.getTargetDim(), T());
        W output_tmp[_f.getTargetDim()];
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          //output[iD] = W();
          output_tmp[iD] = W();
        }

        latticeC[1] = locX;
        latticeC[2] = locY;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (1 - d[1])*(1 - d[2]));
          output_tmp[iD] = W();
        }


        latticeC[1] = locX;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (d[1])*(1 - d[2]));
          output_tmp[iD] = W();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (1 - d[1]) * (1 - d[2]));
          output_tmp[iD] = W();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (d[1]) * (1 - d[2]));
          output_tmp[iD] = W();
        }

        latticeC[1] = locX;
        latticeC[2] = locY;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (1 - d[1]) * (d[2]));
          output_tmp[iD] = W();
        }

        latticeC[1] = locX;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (d[1]) * (d[2]));
          output_tmp[iD] = W();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (1 - d[1]) * (d[2]));
          output_tmp[iD] = W();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (d[1]) * (d[2]));
          output_tmp[iD] = W();
        }
      }
      dataSize += _f.getTargetDim();
      dataFound ++;
      //      for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      //        output[iD] = output_tmp[iD];
      //      }
    }
  }

  if (_communicateToAll) {
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(dataFound, MPI_SUM);
    singleton::mpi().reduceAndBcast(dataSize, MPI_SUM);
#endif
    dataSize /= dataFound;
#ifdef PARALLEL_MODE_MPI
    for (int iD = 0; iD < dataSize; ++iD) {
      singleton::mpi().reduceAndBcast(output[iD], MPI_SUM);
    }
#endif
    for (int iD = 0; iD < dataSize; ++iD) {
      output[iD]/=dataFound;
    }
  } else {
    if (dataFound!=0) {
      dataSize /= dataFound;
      for (int iD = 0; iD < dataSize; ++iD) {
        output[iD]/=dataFound;
      }
    }
  }

  if (dataFound>0) {
    return true;
  }
  return false;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::SuperLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T, T>& f, SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, f.getTargetDim()),
    _f(f)
{
  this->getName() = "fromAnalyticalF";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  this->_sLattice.getCuboidGeometry().getPhysR(physR,input);
  _f(output,physR);
  return true;
}

//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::BlockLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T, T>& f, BlockLattice3D<T, DESCRIPTOR>& sLattice,
  BlockGeometry3D<T>& superGeometry, CuboidGeometry3D<T>& cuboidGeometry)
  : BlockLatticeF3D<T, DESCRIPTOR>(sLattice, f.getTargetDim()),
    _f(f),
    _superGeometry(superGeometry),
    _cuboidGeometry(cuboidGeometry)
{
  this->getName() = "fromAnalyticalF";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  _superGeometry.getPhysR(physR,input[0],input[1],input[2] );
  _f(output,physR);
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
SmoothBlockIndicator3D<T, DESCRIPTOR>::SmoothBlockIndicator3D(
  IndicatorF3D<T>& f, T h )
  : BlockDataF3D<T, T>((int)((f.getMax()[0] - f.getMin()[0] + 4*h) / h + 1.5),
                       (int)((f.getMax()[1] - f.getMin()[1] + 4*h) / h + 1.5), (int)((f.getMax()[2] - f.getMin()[2] + 4*h) / h + 1.5)),
    _f(f),
    _h(h)
{
  this->getName() = "SmoothBlockIndicator3D";

  for (int iX=0; iX<this->getBlockData().getNx(); iX++) {
    for (int iY=0; iY<this->getBlockData().getNy(); iY++) {
      for (int iZ=0; iZ<this->getBlockData().getNz(); iZ++) {
        T value = T();
        bool output[1];
        for (int iPop=0; iPop<DESCRIPTOR<T>::q; iPop++) {
          T input[] = {_f.getMin()[0] - 2*_h + (iX + DESCRIPTOR<T>::c[iPop][0])
                       *_h, _f.getMin()[1] - 2*_h + (iY + DESCRIPTOR<T>::c[iPop][1])*_h,
                       _f.getMin()[2] - 2*_h + (iZ + DESCRIPTOR<T>::c[iPop][2])*_h
                      };
          _f(output,input);
          value += (int)output[0]*DESCRIPTOR<T>::t[iPop];
        }
        //std::cout << iX << "/" << iY << "/" << iZ  << std::endl;
        this->getBlockData().get(iX,iY,iZ,0) = value;
      }
    }
  }
}
/*
template<typename T, template<typename U> class DESCRIPTOR>
bool SmoothBlockIndicator3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  _superGeometry.getPhysR(physR,input[0],input[1],input[2] );
  _f(output,physR);
  return true;
}*/


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
SuperLatticeInterpPhysVelocity3Degree3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& conv, int range)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "Interp3DegreeVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>* foo =
      new BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>(
      sLattice.getExtendedBlockLattice(iC),
      conv,
      &sLattice.getCuboidGeometry().get(this->_sLattice.getLoadBalancer().
                                        glob(iC)),
      sLattice.getOverlap(),
      range);
    bLattices.push_back(foo);
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[], const T input[], const int iC)
{
  bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(output,
      input);
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpPhysVelocity3Degree3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, LBconverter<T>& conv,
  Cuboid3D<T>* c, int overlap, int range)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
    _conv(conv),
    _cuboid(c),
    _overlap(overlap),
    _range(range)
{
  this->getName() = "BlockLatticeInterpVelocity3Degree3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpPhysVelocity3Degree3D(
  const BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>& rhs) :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
  _conv(rhs._conv),
  _cuboid(rhs._cuboid),
  _overlap(rhs._overlap),
  _range(rhs._range)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
void BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[3], const T input[3])
{
  T u[3], rho, volume;
  int latIntPos[3] = {0};
  T latPhysPos[3] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

  latIntPos[0]+=_overlap;
  latIntPos[1] += _overlap;
  latIntPos[2] += _overlap;

  volume=T(1);
  for (int i = -_range; i <= _range+1; ++i) {
    for (int j = -_range; j <= _range+1; ++j) {
      for (int k = -_range; k <= _range+1; ++k) {

        this->_blockLattice.get(latIntPos[0]+i, latIntPos[1]+j,
                                latIntPos[2]+k).computeRhoU(rho, u);
        for (int l = -_range; l <= _range+1; ++l) {
          if (l != i) {
            volume *= (input[0] - (latPhysPos[0]+ l *_cuboid->getDeltaR()))
                      / (latPhysPos[0] + i *_cuboid->getDeltaR()
                         - (latPhysPos[0] + l *_cuboid->getDeltaR()));
          }
        }
        for (int m = -_range; m <= _range+1; ++m) {
          if (m != j) {
            volume *= (input[1]
                       - (latPhysPos[1] + m *_cuboid->getDeltaR()))
                      / (latPhysPos[1] + j * _cuboid->getDeltaR()
                         - (latPhysPos[1] + m * _cuboid->getDeltaR()));
          }
        }
        for (int n = -_range; n <= _range+1; ++n) {
          if (n != k) {
            volume *= (input[2]
                       - (latPhysPos[2] + n * _cuboid->getDeltaR()))
                      / (latPhysPos[2] + k * _cuboid->getDeltaR()
                         - (latPhysPos[2] + n * _cuboid->getDeltaR()));
          }
        }
        output[0] += u[0] * volume;
        output[1] += u[1] * volume;
        output[2] += u[2] * volume;
        volume=T(1);
      }
    }
  }

  output[0] = _conv.physVelocity(output[0]);
  output[1] = _conv.physVelocity(output[1]);
  output[2] = _conv.physVelocity(output[2]);
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::
SuperLatticeInterpDensity3Degree3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& conv, int range)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "Interp3DegreeDensity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>* foo =
      new BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>(
      sLattice.getExtendedBlockLattice(iC),
      conv,
      &sLattice.getCuboidGeometry().get(this->_sLattice.getLoadBalancer().
                                        glob(iC)),
      sLattice.getOverlap(),
      range);
    bLattices.push_back(foo);
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(T output[],
    const T input[], const int iC)
{
  bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(output,
      input);
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpDensity3Degree3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, LBconverter<T>& conv,
  Cuboid3D<T>* c, int overlap, int range)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
    _conv(conv),
    _cuboid(c),
    _overlap(overlap),
    _range(range)
{
  this->getName() = "BlockLatticeInterpDensity3Degree3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::
BlockLatticeInterpDensity3Degree3D(
  const BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>& rhs) :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
  _conv(rhs._conv),
  _cuboid(rhs._cuboid),
  _overlap(rhs._overlap),
  _range(rhs._range)
{
}

template <typename T, template <typename U> class DESCRIPTOR>
void BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(
  T output[DESCRIPTOR<T>::q], const T input[3])
{
  T volume, density;
  int latIntPos[3] = {0};
  T latPhysPos[3] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

  latIntPos[0] += _overlap;
  latIntPos[1] += _overlap;
  latIntPos[2] += _overlap;

  for (unsigned iPop = 0; iPop < DESCRIPTOR<T>::q; ++iPop) {
    volume=T(1);
    for (int i = -_range; i <= _range+1; ++i) {
      for (int j = -_range; j <= _range+1; ++j) {
        for (int k = -_range; k <= _range+1; ++k) {

          density = this->_blockLattice.get(latIntPos[0]+i, latIntPos[1]+j,
                                            latIntPos[2]+k).operator[](iPop);
          for (int l = -_range; l <= _range+1; ++l) {
            if (l != i) {
              volume *= (input[0] - (latPhysPos[0]+ l *_cuboid->getDeltaR()))
                        / (latPhysPos[0] + i *_cuboid->getDeltaR()
                           - (latPhysPos[0] + l *_cuboid->getDeltaR()));
            }
          }
          for (int m = -_range; m <= _range+1; ++m) {
            if (m != j) {
              volume *= (input[1]
                         - (latPhysPos[1] + m *_cuboid->getDeltaR()))
                        / (latPhysPos[1] + j * _cuboid->getDeltaR()
                           - (latPhysPos[1] + m * _cuboid->getDeltaR()));
            }
          }
          for (int n = -_range; n <= _range+1; ++n) {
            if (n != k) {
              volume *= (input[2]
                         - (latPhysPos[2] + n * _cuboid->getDeltaR()))
                        / (latPhysPos[2] + k * _cuboid->getDeltaR()
                           - (latPhysPos[2] + n * _cuboid->getDeltaR()));
            }
          }
          output[iPop] += density * volume;

          volume=T(1);
        }
      }
    }
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::SuperLatticeSmoothDiracDelta3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& sGeometry,
  LBconverter<T>& conv) :
  SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "SuperLatticeSmoothDiracDelta3D";
  int maxC = this->_sLattice.getCuboidGeometry().getNc();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {

    if (sLattice.getLoadBalancer().rank(iC) == singleton::mpi().getRank()) {
      int locIC = sLattice.getLoadBalancer().loc(iC);

      BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>* foo =
        new BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>(
        sLattice.getExtendedBlockLattice(locIC),
        sGeometry.getBlockGeometry(locIC),
        &sLattice.getCuboidGeometry().get(iC),
        conv, sLattice.getOverlap());
      bLattices.push_back(foo);
    }
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::operator()(T delta[4][4][4],
    const T physPos[3], const int iC)
{
  bLattices[this->_sLattice.getLoadBalancer().loc(iC)]->operator()(delta,
      physPos);
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::BlockLatticeSmoothDiracDelta3D(
  BlockLattice3D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>&
  blockGeometry, Cuboid3D<T>* cuboid, LBconverter<T>& conv, int overlap)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3), _cuboid(cuboid),
    _conv(conv),  _bGeometry(blockGeometry), _overlap(overlap)
{
  this->getName() = "BlockLatticeSmoothDiracDelta3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::BlockLatticeSmoothDiracDelta3D(
  const BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>& rhs)
  :
  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3), _cuboid(rhs._cuboid),
  _conv(rhs._conv), _bGeometry(rhs._bGeometry), _overlap(rhs._overlap)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
void BlockLatticeSmoothDiracDelta3D<T, DESCRIPTOR>::operator()(
  T delta[4][4][4], const T physPos[])
{
  int range = 1;
  T a, b, c = T();
  int latticeRoundedPosP[3] = { 0 };
  T physRoundedPosP[3] = { T() };
  T physLatticeL = _conv.getLatticeL();

  T counter = 0.;

  _cuboid->getLatticeR(latticeRoundedPosP, physPos);
  _cuboid->getPhysR(physRoundedPosP, latticeRoundedPosP);

  for (int i = -range; i <= range + 1; ++i) {
    for (int j = -range; j <= range + 1; ++j) {
      for (int k = -range; k <= range + 1; ++k) {
        delta[i+range][j+range][k+range] = T(1);
        // a, b, c in lattice units cause physical ones get cancelled
        a = (physRoundedPosP[0] + i * physLatticeL - physPos[0])
            / physLatticeL;
        b =  (physRoundedPosP[1] + j * physLatticeL - physPos[1])
             / physLatticeL;
        c = (physRoundedPosP[2] + k * physLatticeL - physPos[2])
            / physLatticeL;

        // the for loops already define that a, b, c are smaller than 2
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + cos(M_PI * a / 2.));
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + cos(M_PI * b / 2.));
        delta[i+range][j+range][k+range] *= 1. / 4 * (1 + cos(M_PI * c / 2.));

        counter += delta[i+range][j+range][k+range];
      }
    }
  }

  //  if (!util::nearZero(counter - T(1))){
  //    // sum of delta has to be one
  //    std::cout << "[" << this->getName() << "] " <<
  //        "Delta summed up does not equal 1 but = " <<
  //        counter << std::endl;
  //  }

}

}  // end namespace olb

#endif
