/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef INTERPOLATION_F_2D_HH
#define INTERPOLATION_F_2D_HH


#include "functors/interpolationF2D.h"
#include "core/superLattice2D.h"
#include "dynamics/lbHelpers.h"


namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
AnalyticalFfromSuperLatticeF2D<T,DESCRIPTOR>::
AnalyticalFfromSuperLatticeF2D(SuperLatticeF2D<T,DESCRIPTOR>& f, bool communicateToAll,
                               int overlap)
  : AnalyticalF2D<T,T>(f.getTargetDim()), _f(f),
    _cg(_f.getSuperLattice().getCuboidGeometry()), _communicateToAll(communicateToAll),
    _overlap(overlap)
{
  this->getName() = "fromSuperLatticeF";
  if (overlap==-1) {
    _overlap = f.getSuperLattice().getOverlap();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
bool AnalyticalFfromSuperLatticeF2D<T,DESCRIPTOR>::operator() (T output[], const T physC[])
{
  int latticeR[3] = {0,0,0};
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = T();
  }
  if (!(_cg.getLatticeR(latticeR,physC))) {
    return false;
  }

  // convert to lattice coordinates
  T d[2];

  _f.getSuperLattice().communicate();

  int locX, locY;

  int dataSize = 0;
  int dataFound = 0;

  int latticeC[3]= {};
  T physRiC[2];


  for (int iC = 0; iC < _f.getSuperLattice().getLoadBalancer().size(); ++iC) {
    latticeC[0] = _f.getSuperLattice().getLoadBalancer().glob(iC);
    _cg.get(latticeC[0]).getFloorLatticeR(latticeR, physC);
    if (latticeR[0] >= -_overlap && latticeR[0] + 1 < _cg.get(latticeC[0]).getNx() + _overlap &&
        latticeR[1] >= -_overlap && latticeR[1] + 1 < _cg.get(latticeC[0]).getNy() + _overlap ) {

      locX = latticeR[0];
      locY = latticeR[1];

      _cg.get(latticeC[0]).getPhysR(physRiC, locX,locY);

      d[0] = (physC[0] - physRiC[0]);
      d[1] = (physC[1] - physRiC[1]);

      d[0] /= _cg.get(latticeC[0]).getDeltaR();
      d[1] /= _cg.get(latticeC[0]).getDeltaR();

      T outputTmp[_f.getTargetDim()];
      for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
        //output[iD] = T();
        outputTmp[iD] = T();
      }

      latticeC[1] = locX;
      latticeC[2] = locY;
      _f(outputTmp,latticeC);
      for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
        output[iD] += (outputTmp[iD]*(1-d[0])*(1-d[1]));
        outputTmp[iD] = T();
      }
      latticeC[1] = locX;
      latticeC[2] = locY+1;
      _f(outputTmp,latticeC);
      for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
        output[iD] += (outputTmp[iD]*(1-d[0])*(d[1]));
        outputTmp[iD] = T();
      }
      latticeC[1] = locX+1;
      latticeC[2] = locY;
      _f(outputTmp,latticeC);
      for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
        output[iD] += (outputTmp[iD]*(d[0])*(1-d[1]));
        outputTmp[iD] = T();
      }
      latticeC[1] = locX+1;
      latticeC[2] = locY+1;
      _f(outputTmp,latticeC);
      for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
        output[iD] += (outputTmp[iD]*(d[0])*(d[1]));
        outputTmp[iD] = T();
      }

      dataSize += _f.getTargetDim();
      dataFound ++;
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



template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR>::
SuperLatticeFfromAnalyticalF2D(AnalyticalF2D<T,T>& f, SuperLattice2D<T,DESCRIPTOR>& sLattice,
                               SuperGeometry2D<T>& sg)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, f.getTargetDim()), _f(f), _sg(sg)
{
  this->getName() = "fromAnalyticalF";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  // convert to physical coordinates
  T physCoordinate[2] = {};
  physCoordinate[0] = _sg.getPhysR(input[0],input[1],input[2])[0];
  physCoordinate[1] = _sg.getPhysR(input[0],input[1],input[2])[1];

  _f(output,physCoordinate);

//std::cout << physCoordinate[0] << "/"<<physCoordinate[1] << "/"<<output[0]<< "/"<<output[1] <<std::endl;

  return true;
}



} // end namespace olb

#endif
