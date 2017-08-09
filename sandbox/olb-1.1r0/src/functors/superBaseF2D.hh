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

#ifndef SUPER_BASE_F_2D_HH
#define SUPER_BASE_F_2D_HH


#include "superBaseF2D.h"

namespace olb {

template<typename T, typename W>
SuperF2D<T,W>::SuperF2D(SuperStructure2D<T>& superStructure, int targetDim)
  : GenericF<W,int>(targetDim,3), _superStructure(superStructure) { }

template<typename T, typename W>
SuperF2D<T,W>::~SuperF2D()
{
  for (auto&& element : _blockF) {
    delete element;
  }
  _blockF.resize(0);
}

template<typename T, typename W>
SuperStructure2D<T>& SuperF2D<T,W>::getSuperStructure()
{
  return _superStructure;
}

template<typename T, typename W>
BlockF2D<T>& SuperF2D<T,W>::getBlockF(int iCloc)
{
  return *(_blockF[iCloc]);
}

template <typename T,typename BaseType>
SuperDataF2D<T,BaseType>::SuperDataF2D(SuperData2D<T,BaseType>& superData)
  : SuperF2D<T,BaseType>( superData, superData.getDataSize() ), _superData(superData)
{}

template <typename T,typename BaseType>
bool SuperDataF2D<T,BaseType>::operator() (T output[], const int input[])
{
  for (int i=0; i < _superData.getDataSize(); i++) {
    output[i] = _superData.get(_superData.getLoadBalancer().loc(input[0]))
                .get(input[1]+_superData.getOverlap(), input[2]+_superData.getOverlap(), i);
  }
  return true;
}

template <typename T,typename BaseType>
SuperData2D<T,BaseType>& SuperDataF2D<T,BaseType>::getSuperData()
{
  return _superData;
}


template <typename T, typename W>
SuperIdentity2D<T,W>::SuperIdentity2D(SuperF2D<T,W>& f)
  : SuperF2D<T,W>(f.getSuperStructure() ,f.getTargetDim() ), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename W>
bool SuperIdentity2D<T,W>::operator()(W output[], const int input[])
{
  _f(output, input);
  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF2D<T,DESCRIPTOR>::SuperLatticeF2D(SuperLattice2D<T,DESCRIPTOR>& superLattice,
    int targetDim)
  : SuperF2D<T,T>(superLattice, targetDim), _sLattice(superLattice) { }

template <typename T, template <typename U> class DESCRIPTOR>
SuperLattice2D<T,DESCRIPTOR>& SuperLatticeF2D<T,DESCRIPTOR>::getSuperLattice()
{
  return _sLattice;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysF2D<T,DESCRIPTOR>::SuperLatticePhysF2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter,
 int targetDim)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, targetDim), _converter(converter) { }

template <typename T, template <typename U> class DESCRIPTOR>
LBconverter<T> const& SuperLatticePhysF2D<T,DESCRIPTOR>::getConverter() const
{
  return this->_converter;
}


} // end namespace olb

#endif
