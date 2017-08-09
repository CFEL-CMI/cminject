/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2008 Jonas Latt
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
 * Dynamics for a generic 2D block lattice view -- generic implementation.
 */
#ifndef BLOCK_LATTICE_VIEW_2D_HH
#define BLOCK_LATTICE_VIEW_2D_HH

#include "blockLatticeView2D.h"
#include "cell.h"

namespace olb {

////////////////////// Class BlockLatticeView2D /////////////////////////

template<typename T, template<typename U> class Lattice>
BlockLatticeView2D<T,Lattice>::BlockLatticeView2D (
  BlockLatticeStructure2D<T,Lattice>& originalLattice_ )
  : BlockLatticeStructure2D<T,Lattice>(originalLattice_.getNx(),originalLattice_.getNy()), originalLattice(&originalLattice_),
    x0(0), y0(0)
{ }

template<typename T, template<typename U> class Lattice>
BlockLatticeView2D<T,Lattice>::BlockLatticeView2D (
  BlockLatticeStructure2D<T,Lattice>& originalLattice_,
  int x0_, int x1_, int y0_, int y1_ )
  : BlockLatticeStructure2D<T,Lattice>(x1_-x0_+1,y1_-y0_+1), originalLattice(&originalLattice_),
    x0(x0_), y0(y0_)
{
  OLB_PRECONDITION(x1_ < originalLattice->getNx());
  OLB_PRECONDITION(x0_ <= x1_);
  OLB_PRECONDITION(y1_ < originalLattice->getNy());
  OLB_PRECONDITION(y0_ <= y1_);
}

template<typename T, template<typename U> class Lattice>
BlockLatticeView2D<T,Lattice>::~BlockLatticeView2D()
{
}


template<typename T, template<typename U> class Lattice>
BlockLatticeView2D<T,Lattice>::BlockLatticeView2D(BlockLatticeView2D<T,Lattice> const& rhs)
  : BlockLatticeStructure2D<T,Lattice>(rhs._nx,rhs._ny), originalLattice(rhs.originalLattice),
    x0(rhs.x0), y0(rhs.y0)
{ }

template<typename T, template<typename U> class Lattice>
BlockLatticeView2D<T,Lattice>& BlockLatticeView2D<T,Lattice>::operator= (
  BlockLatticeView2D<T,Lattice> const& rhs )
{
  BlockLatticeView2D<T,Lattice> tmp(rhs);
  swap(tmp);
  return *this;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::swap (
  BlockLatticeView2D<T,Lattice>& rhs)
{
  std::swap(this->_nx, rhs._nx);
  std::swap(this->_ny, rhs._ny);
  std::swap(x0, rhs.x0);
  std::swap(y0, rhs.y0);
  std::swap(originalLattice, rhs.originalLattice);
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& BlockLatticeView2D<T,Lattice>::get(int iX, int iY)
{
  OLB_PRECONDITION(iX<this->_nx);
  OLB_PRECONDITION(iY<this->_ny);
  return originalLattice->get(iX+x0, iY+y0);
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& BlockLatticeView2D<T,Lattice>::get (
  int iX, int iY ) const
{
  OLB_PRECONDITION(iX<this->_nx);
  OLB_PRECONDITION(iY<this->_ny);
  return originalLattice->get(iX+x0, iY+y0);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::initialize()
{
  originalLattice->initialize();
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::defineDynamics (
  int x0_, int x1_, int y0_, int y1_,
  Dynamics<T,Lattice>* dynamics )
{
  originalLattice->defineDynamics( x0_+x0, x1_+x0,
                                   y0_+y0, y1_+y0,
                                   dynamics );

}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::defineDynamics (
  int iX, int iY, Dynamics<T,Lattice>* dynamics )
{
  originalLattice->defineDynamics( iX+x0, iY+y0, dynamics );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::specifyStatisticsStatus (
  int x0_, int x1_, int y0_, int y1_, bool status )
{
  originalLattice->specifyStatisticsStatus(
    x0_+x0, x1_+x0,
    y0_+y0, y1_+y0,
    status );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::collide (
  int x0_, int x1_, int y0_, int y1_ )
{
  originalLattice->collide( x0_+x0, x1_+x0,
                            y0_+y0, y1_+y0 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::collide()
{
  originalLattice->collide( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}
/*
template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::staticCollide (
  int x0_, int x1_, int y0_, int y1_,
  TensorFieldBase2D<T,2> const& u)
{
  originalLattice->staticCollide( x0_+x0, x1_+x0,
                                  y0_+y0, y1_+y0, u );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::staticCollide (
  TensorFieldBase2D<T,2> const& u )
{
  originalLattice->staticCollide( x0, x0+this->_nx-1, y0, y0+this->_ny-1, u);
}
*/

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::stream(int x0_, int x1_, int y0_, int y1_)
{
  originalLattice->stream(x0_+x0, x1_+x0, y0_+y0, y1_+y0);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::stream(bool periodic)
{
  OLB_PRECONDITION(!periodic);
  originalLattice->stream( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
  postProcess();
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::collideAndStream(int x0_, int x1_, int y0_, int y1_)
{
  originalLattice->collideAndStream(x0_+x0, x1_+x0, y0_+y0, y1_+y0);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::collideAndStream(bool periodic)
{
  OLB_PRECONDITION(!periodic);
  originalLattice->collideAndStream( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
  postProcess();
}

template<typename T, template<typename U> class Lattice>
T BlockLatticeView2D<T,Lattice>::computeAverageDensity (
  int x0_, int x1_, int y0_, int y1_ ) const
{
  return originalLattice->computeAverageDensity( x0_+x0, x1_+x0,
         y0_+y0, y1_+y0 );
}

template<typename T, template<typename U> class Lattice>
T BlockLatticeView2D<T,Lattice>::computeAverageDensity() const
{
  return originalLattice->computeAverageDensity( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::stripeOffDensityOffset (
  int x0_, int x1_, int y0_, int y1_, T offset )
{
  originalLattice->stripeOffDensityOffset( x0_+x0, x1_+x0,
      y0_+y0, y1_+y0, offset );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::stripeOffDensityOffset(T offset)
{
  originalLattice->stripeOffDensityOffset(x0, x0+this->_nx-1, y0, y0+this->_ny-1, offset);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::forAll (
  int x0_, int x1_, int y0_, int y1_, WriteCellFunctional<T,Lattice> const& application )
{
  originalLattice->forAll( x0_+x0, x1_+x0, y0_+y0, y1_+y0, application );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::forAll(WriteCellFunctional<T,Lattice> const& application)
{
  originalLattice->forAll(x0, x0+this->_nx-1, y0, y0+this->_ny-1, application);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::addPostProcessor (
  PostProcessorGenerator2D<T,Lattice> const& ppGen)
{
  PostProcessorGenerator2D<T,Lattice>* shiftedGen = ppGen.clone();
  shiftedGen->shift(x0,y0);
  originalLattice->addPostProcessor(*shiftedGen);
  delete shiftedGen;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::resetPostProcessors()
{
  originalLattice->resetPostProcessors();
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::postProcess()
{
  originalLattice -> postProcess(x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::postProcess (
  int x0_, int x1_, int y0_, int y1_ )
{
  originalLattice -> postProcess( x0_+x0, x1_+x0, y0_+y0, y1_+y0 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::addLatticeCoupling (
  LatticeCouplingGenerator2D<T,Lattice> const& lcGen,
  std::vector<SpatiallyExtendedObject2D*> partners )
{
  LatticeCouplingGenerator2D<T,Lattice>* shiftedGen = lcGen.clone();
  shiftedGen->shift(x0,y0);
  originalLattice->addLatticeCoupling(*shiftedGen, partners);
  delete shiftedGen;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::executeCoupling()
{
  originalLattice -> executeCoupling(x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::executeCoupling (
  int x0_, int x1_, int y0_, int y1_ )
{
  originalLattice -> executeCoupling( x0_+x0, x1_+x0, y0_+y0, y1_+y0 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView2D<T,Lattice>::
subscribeReductions(Reductor<T>& reductor)
{
  originalLattice -> subscribeReductions(reductor);
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T>& BlockLatticeView2D<T,Lattice>::getStatistics()
{
  return originalLattice->getStatistics();
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T> const& BlockLatticeView2D<T,Lattice>::getStatistics() const
{
  return originalLattice->getStatistics();
}

template<typename T, template<typename U> class Lattice>
SpatiallyExtendedObject2D* BlockLatticeView2D<T,Lattice>::getComponent(int iBlock)
{
  OLB_PRECONDITION( iBlock==0 );
  return this;
}

template<typename T, template<typename U> class Lattice>
SpatiallyExtendedObject2D const* BlockLatticeView2D<T,Lattice>::getComponent(int iBlock) const
{
  OLB_PRECONDITION( iBlock==0 );
  return this;
}

template<typename T, template<typename U> class Lattice>
multiPhysics::MultiPhysicsId BlockLatticeView2D<T,Lattice>::getMultiPhysicsId() const
{
  return multiPhysics::getMultiPhysicsBlockId<T,Lattice>();
}

}  // namespace olb

#endif
