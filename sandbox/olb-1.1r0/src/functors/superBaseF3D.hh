/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2016 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *                          Albert Mink, Benjamin Förster
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

#ifndef SUPER_BASE_F_3D_HH
#define SUPER_BASE_F_3D_HH


#include "superBaseF3D.h"

namespace olb {

template <typename T, typename W>
SuperF3D<T,W>::SuperF3D(SuperStructure3D<T>& superStructure, int targetDim)
  : GenericF<W,int>(targetDim,4), _superStructure(superStructure) { }

template <typename T, typename W>
SuperF3D<T,W>::~SuperF3D()
{
  for (auto&& element : _blockF) {
    delete element;
  }
  _blockF.resize(0);
}

template <typename T, typename W>
SuperStructure3D<T>& SuperF3D<T,W>::getSuperStructure()
{
  return _superStructure;
}

template <typename T, typename W>
BlockF3D<T>& SuperF3D<T,W>::getBlockF(int iCloc)
{
  return *(_blockF[iCloc]);
}




template <typename T,typename BaseType>
SuperDataF3D<T,BaseType>::SuperDataF3D(SuperData3D<T,BaseType>& superData)
  : SuperF3D<T,BaseType>( superData, superData.getDataSize() ), _superData(superData)
{
}


template <typename T,typename BaseType>
bool SuperDataF3D<T,BaseType>::operator() (BaseType output[], const int input[])
{
  int iC_loc = _superData.getLoadBalancer().loc(input[0]);
  for (int i=0; i < _superData.getDataSize(); i++) {
    output[i] = _superData.get(iC_loc).get(input[1] + _superData.getOverlap(),
                                           input[2] + _superData.getOverlap(),
                                           input[3] + _superData.getOverlap(), i);
  }
  return true;
}

template <typename T,typename BaseType>
SuperData3D<T,BaseType>& SuperDataF3D<T,BaseType>::getSuperData()
{
  return _superData;
}



template <typename T, typename W>
SuperIdentity3D<T,W>::SuperIdentity3D(SuperF3D<T,W>& f)
  : SuperF3D<T,W>(f.getSuperStructure() ,f.getTargetDim() ), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename W>
bool SuperIdentity3D<T,W>::operator()(W output[], const int input[])
{
  _f(output, input);
  return true;
}



template <typename T, typename W>
SuperIdentityOnSuperIndicatorF3D<T,W>::SuperIdentityOnSuperIndicatorF3D(SuperF3D<T,W>& f,
    SuperIndicatorF3D<T>& indicatorF,
    W defaultValue)
  : SuperF3D<T,W>(f.getSuperStructure() ,f.getTargetDim() ),
    _f(f),
    _indicatorF(indicatorF),
    _defaultValue(defaultValue)
{
  this->getName() = _f.getName() + std::string("_on_") + _indicatorF.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename W>
bool SuperIdentityOnSuperIndicatorF3D<T,W>::operator()(W output[], const int input[])
{
  bool indic;
  _indicatorF(&indic, input);
  if (indic) {
    _f(output, input);
  } else {
    for (int i=0; i<_f.getTargetDim(); i++) {
      output[i] = _defaultValue;
    }
  }
  return true;
}


template <typename T, template <typename U> class Lattice>
SuperLatticeF3D<T,Lattice>::SuperLatticeF3D(SuperLattice3D<T,Lattice>& superLattice,
    int targetDim)
  : SuperF3D<T,T>(superLattice, targetDim), _sLattice(superLattice) { }

template <typename T, template <typename U> class Lattice>
SuperLattice3D<T,Lattice>& SuperLatticeF3D<T,Lattice>::getSuperLattice()
{
  return _sLattice;
}

template <typename T, template <typename U> class Lattice>
SuperLatticePhysF3D<T,Lattice>::SuperLatticePhysF3D
(SuperLattice3D<T,Lattice>& sLattice, const LBconverter<T>& converter,
 int targetDim)
  : SuperLatticeF3D<T,Lattice>(sLattice, targetDim), _converter(converter) { }

template <typename T, template <typename U> class Lattice>
LBconverter<T> const& SuperLatticePhysF3D<T,Lattice>::getConverter() const
{
  return this->_converter;
}


template <typename T, template <typename U> class Lattice>
ComposedSuperLatticeF3D<T,Lattice>::ComposedSuperLatticeF3D
(SuperLatticeF3D<T,Lattice>& f0, SuperLatticeF3D<T,Lattice>& f1,
 SuperLatticeF3D<T,Lattice>& f2)
  : SuperLatticeF3D<T,Lattice>(f0.getSuperLattice(), 3), _f0(f0), _f1(f1), _f2(f2)
{
  this->getName() = "composedSuperLatticeF3D";
}

template <typename T, template <typename U> class Lattice>
bool ComposedSuperLatticeF3D<T,Lattice>::operator() (T output[], const int input[])
{
  T tmp[3] = {};
  SuperIdentity3D<T,T> ff0(_f0), ff1(_f1), ff2(_f2);
  _f0(tmp,input);
  output[0]=tmp[0];
  _f1(tmp,input);
  output[1]=tmp[0];
  _f2(tmp,input);
  output[2]=tmp[0];
  return true;
}


} // end namespace olb

#endif
