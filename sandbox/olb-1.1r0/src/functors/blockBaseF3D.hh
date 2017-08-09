/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_BASE_F_3D_HH
#define BLOCK_BASE_F_3D_HH


#include "functors/blockBaseF3D.h"


namespace olb {


template <typename T>
BlockF3D<T>::BlockF3D(BlockStructure3D& blockStructure, int targetDim)
  : GenericF<T,int>(targetDim,3), _blockStructure(blockStructure) { }

template <typename T>
BlockStructure3D& BlockF3D<T>::getBlockStructure() const
{
  return _blockStructure;
}

template <typename T>
std::vector<T> BlockF3D<T>::getMinValue()
{
  T min[this->getTargetDim()];
  T minTmp[this->getTargetDim()];
  this->operator()(min,0,0);
  for (int iX = 1; iX < _blockStructure.getNx(); ++iX) {
    for (int iY = 1; iY < _blockStructure.getNy(); ++iY) {
      for (int iZ = 1; iZ < _blockStructure.getNz(); ++iZ) {
        this->operator()(minTmp,iX,iY,iZ);
        for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
          if (min[iDim] > minTmp[iDim] ) {
            min[iDim] = minTmp[iDim];
          }
        }
      }
    }
  }
  std::vector<T> minV(min,min+this->getTargetDim());
  return minV;
}


template <typename T>
std::vector<T> BlockF3D<T>::getMaxValue()
{
  T max[this->getTargetDim()];
  T maxTmp[this->getTargetDim()];
  this->operator()(max,0,0);
  for (int iX = 1; iX < _blockStructure.getNx(); ++iX) {
    for (int iY = 1; iY < _blockStructure.getNy(); ++iY) {
      for (int iZ = 1; iZ < _blockStructure.getNz(); ++iZ) {
        this->operator()(maxTmp,iX,iY,iZ);
        for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
          if (max[iDim] > maxTmp[iDim] ) {
            max[iDim] = maxTmp[iDim];
          }
        }
      }
    }
  }
  std::vector<T> maxV(max,max+this->getTargetDim());
  return maxV;
}

template <typename T>
BlockIdentity3D<T>::BlockIdentity3D(BlockF3D<T>& f)
  : BlockF3D<T>(f.getBlockStructure() ,f.getTargetDim() ), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

/*
template <typename T>
BlockIdentity3D<T> BlockIdentity3D<T>::operator=(BlockF3D<T>& rhs)
{
  BlockIdentity3D<T> tmp(rhs);
  return tmp;
}
*/


template <typename T>
bool BlockIdentity3D<T>::operator()(T output[], const int input[])
{
  _f(output,input);
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>::BlockLatticeF3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockStructure, int targetDim)
  : BlockF3D<T>(blockStructure, targetDim), _blockLattice(blockStructure)
{ }
/*
template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>::BlockLatticeF3D(BlockLatticeF3D<T,DESCRIPTOR> const& rhs)
  : BlockF3D<T>(rhs.getBlockStructure(), rhs.getTargetDim() ), _blockLattice(rhs.getBlockLattice())
{ }

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>& BlockLatticeF3D<T,DESCRIPTOR>::operator=(BlockLatticeF3D<T,DESCRIPTOR> const& rhs)
{
  BlockLatticeF3D<T,DESCRIPTOR> tmp(rhs);
  return tmp;
}
*/
template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeStructure3D<T,DESCRIPTOR>& BlockLatticeF3D<T, DESCRIPTOR>::getBlockLattice()
{
  return _blockLattice;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysF3D<T,DESCRIPTOR>::BlockLatticePhysF3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter,
 int targetDim)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, targetDim), _converter(converter) { }






} // end namespace olb

#endif
