/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2015 Mathias J. Krause, Jonas Latt, Patrick Nathen
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
 * BGK Dynamics with adjusted omega -- generic implementation.
 */
#ifndef SMAGORINSKY_BGK_DYNAMICS_HH
#define SMAGORINSKY_BGK_DYNAMICS_HH

#include "smagorinskyBGKdynamics.h"
#include "core/cell.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "core/units.h"
#include "core/units.hh"
#include "math.h"

namespace olb {

///////////////////////// ADM BGK /////////////////////////////

template<typename T, template<typename U> class Lattice>
ADMBGKdynamics<T,Lattice>::ADMBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
void ADMBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics )
{
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *cell.getExternal(rhoIsAt), cell.getExternal(velocityBeginsAt), omega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(*cell.getExternal(rhoIsAt), uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void ADMBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *cell.getExternal(rhoIsAt), u, omega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(*cell.getExternal(rhoIsAt), uSqr);
  }
}


///////////////////////// ForcedADM BGK /////////////////////////////

template<typename T, template<typename U> class Lattice>
ForcedADMBGKdynamics<T,Lattice>::ForcedADMBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
void ForcedADMBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  OstreamManager clout(std::cout,"Forced ADM collide:");
  T rho, u[Lattice<T>::d], utst[Lattice<T>::d];

// this->momenta.computeAllMomenta(cell, rho, utst, pi);

  T* rho_fil = cell.getExternal(filRhoIsAt);
  T* u_filX = cell.getExternal(localFilVelXBeginsAt);
  T* u_filY = cell.getExternal(localFilVelYBeginsAt);
  T* u_filZ = cell.getExternal(localFilVelZBeginsAt);

  u[0] = *u_filX;/// *rho_fil;
  u[1] = *u_filY;/// *rho_fil;
  u[2] = *u_filZ;/// *rho_fil;

  T* force = cell.getExternal(forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *rho_fil, u, omega);

  lbHelpers<T,Lattice>::addExternalForce(cell, u, omega);

  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void ForcedADMBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d];


  T* rho_fil = cell.getExternal(filRhoIsAt);
  T* u_filX = cell.getExternal(localFilVelXBeginsAt);
  T* u_filY = cell.getExternal(localFilVelYBeginsAt);
  T* u_filZ = cell.getExternal(localFilVelZBeginsAt);

  uTemp[0] = *u_filX;/// *rho_fil;
  uTemp[1] = *u_filY;/// *rho_fil;
  uTemp[2] = *u_filZ;/// *rho_fil;


  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *rho_fil, uTemp, omega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}



///////////////////// Class ConSmagorinskyBGKdynamics //////////////////////////
//Consistent Smagorinsky BGK --> Malaspinas/Sagaut

template<typename T, template<typename U> class Lattice>
ConSmagorinskyBGKdynamics<T,Lattice>::ConSmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_)), dx(dx_), dt(dt_)
{ }

template<typename T, template<typename U> class Lattice>
void ConSmagorinskyBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics)
{
  T rho;
  T u[Lattice<T>::d];
  T pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T omega = this->getOmega();

  /**************************************************************/
  T H[util::TensorVal<Lattice<T> >::n];
  T conSmagoR[Lattice<T>::q];
  T S[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T tau_mol = 1./omega;
  T cs2 = 1./Lattice<T>::invCs2;

  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(2*PiNeqNormSqr);

  //Strain rate Tensor
  if (PiNeqNorm != 0) {
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] =
        //Pi
        //-1/(2*cs2*rho)*pi[n];
        //Dx
        //(-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst*smagoConst*dx*dx)*rho*PiNeqNorm))/(smagoConst*smagoConst*dx*dx*rho*PiNeqNorm))*pi[n];
        //Malas
        (-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst*smagoConst)*rho*PiNeqNorm))/(smagoConst*smagoConst*rho*PiNeqNorm))*pi[n];
      //all
      //(-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst*smagoConst*dx*dx*dt)*rho*PiNeqNorm))/(smagoConst*smagoConst*dx*dx*dt*rho*PiNeqNorm))*pi[n];
    }
  } else {
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] = 0;
    }
  }

  //Strain rate Tensor Norm
  T SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
  }
  T SNorm    = sqrt(2*SNormSqr);

  //consistent Samagorinsky additional R term
  for (int q = 0; q < Lattice<T>::q; ++q) {
    T t = Lattice<T>::t[q]; //lattice weights

    //Hermite-Polynom H = c*c-cs^2*kronDelta
    H[0] = Lattice<T>::c[q][0]*Lattice<T>::c[q][0]-cs2;
    H[1] = Lattice<T>::c[q][0]*Lattice<T>::c[q][1];
    H[2] = Lattice<T>::c[q][1]*Lattice<T>::c[q][1]-cs2;//2D
    if (util::TensorVal<Lattice<T> >::n == 6) {
      H[2] = Lattice<T>::c[q][0]*Lattice<T>::c[q][2];//3D
      H[3] = Lattice<T>::c[q][1]*Lattice<T>::c[q][1]-cs2;
      H[4] = Lattice<T>::c[q][1]*Lattice<T>::c[q][2];
      H[5] = Lattice<T>::c[q][2]*Lattice<T>::c[q][2]-cs2;
    }

    //contraction or scalar product H*S
    T contractHS = H[0]*S[0] + 2.0*H[1]*S[1] + H[2]*S[2];
    if (util::TensorVal<Lattice<T> >::n == 6) {
      contractHS += H[2]*S[2] + H[3]*S[3] + 2.0*H[4]*S[4] + H[5]*S[5];
    }

    //additional term
    conSmagoR[q] = t*preFactor*SNorm*contractHS;
  }
  /**************************************************************/

  // T uSqr = lbHelpers<T,Lattice>::bgkCollision2(cell, rho, u, omega, conSmagoR);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, omega);
  for (int q = 0; q < Lattice<T>::q; ++q) {
    cell[q] += conSmagoR[q];
  }
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
T ConSmagorinskyBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  //Dx
  //return (T) 0.5*(smagoConst*smagoConst*dx*dx)*Lattice<T>::invCs2*Lattice<T>::invCs2*omega;
  //Malas
  return (T)smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
  //all
  //return (T) 0.5*(smagoConst*smagoConst*dx*dx)*Lattice<T>::invCs2*Lattice<T>::invCs2*dt*omega;
}

template<typename T, template<typename U> class Lattice>
T ConSmagorinskyBGKdynamics<T,Lattice>::computeOmega (T omega0, T preFactor_,
    T rho, T pi[util::TensorVal<Lattice<T> >::n] )
{
  return 0;
}

template<typename T, template<typename U> class Lattice>
void ConSmagorinskyBGKdynamics<T,Lattice>::staticCollide(Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void ConSmagorinskyBGKdynamics<T,Lattice>::setOmega(T omega_)
{
//  this->omega = omega_;
  BGKdynamics<T,Lattice>::setOmega(omega_);
  preFactor = computePreFactor(omega_, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T ConSmagorinskyBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}



////////////////////// Class ConStrainSmagorinskyBGKdynamics //////////////////////////
//Consistent Strain BGK --> Malaspinas/Sagaut
template<typename T, template<typename U> class Lattice>
ConStrainSmagorinskyBGKdynamics<T,Lattice>::ConStrainSmagorinskyBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void ConStrainSmagorinskyBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics)
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
T ConStrainSmagorinskyBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}

template<typename T, template<typename U> class Lattice>
T ConStrainSmagorinskyBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor_,
    T rho, T pi[util::TensorVal<Lattice<T> >::n])
{
  T S[util::TensorVal<Lattice<T> >::n];
  T cs2 = 1./Lattice<T>::invCs2;
  T tau_mol = 1./omega0;

  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.0*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(2*PiNeqNormSqr);

  //Strain Tensor
  if ( !util::nearZero(PiNeqNorm) ) {
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] =
        (-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst*smagoConst)*rho*PiNeqNorm))/(smagoConst*smagoConst*rho*PiNeqNorm))*pi[n];
    }
  } else {
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] = 0;
    }
  }

  //Strain Tensor Norm
  T SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
  }
  T SNorm    = sqrt(2*SNormSqr);

  /// Turbulent realaxation time
  T tau_turb = preFactor_*SNorm;
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

template<typename T, template<typename U> class Lattice>
void ConStrainSmagorinskyBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void ConStrainSmagorinskyBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  BGKdynamics<T,Lattice>::setOmega(omega_);
  preFactor = computePreFactor(omega_, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T ConStrainSmagorinskyBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

///////////////////////// DYNAMIC SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
DynSmagorinskyBGKdynamics<T,Lattice>::DynSmagorinskyBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_, T dx_, T dt_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), dx(dx_), dt(dt_)
{ }

template<typename T, template<typename U> class Lattice>
void DynSmagorinskyBGKdynamics<T,Lattice>::collide( Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T preFactor_ = computePreFactor(this->getOmega(),
                                  *cell.getExternal(Lattice<T>::ExternalField::smagoConstIsAt));
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor_, rho, pi, cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void DynSmagorinskyBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T preFactor_ = computePreFactor(this->getOmega(),
                                  *cell.getExternal(Lattice<T>::ExternalField::smagoConstIsAt));
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor_, rho, pi, cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void DynSmagorinskyBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  BGKdynamics<T,Lattice>::setOmega(omega_);
// this->omega = omega_;
}

template<typename T, template<typename U> class Lattice>
T DynSmagorinskyBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T preFactor_ = computePreFactor(this->getOmega(),
                                  *cell.getExternal(Lattice<T>::ExternalField::smagoConstIsAt));
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor_, rho, pi, cell);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T DynSmagorinskyBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}

template<typename T, template<typename U> class Lattice>
T DynSmagorinskyBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor_, T rho,
    T pi[util::TensorVal<Lattice<T> >::n],Cell<T,Lattice>& cell )
{
  // computation of the relaxation time
  T v_t = 0;
  T* dynSmago = cell.getExternal(smagoConstIsAt);
  //clout << "dynSmago ... " << *dynSmago<<std::endl;
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  //cout << "dynSmago"<< dynSmago<endl;
  v_t = *dynSmago*dx*dx*PiNeqNorm;
  T tau_t = 3*v_t;
  T tau_0 = 1/omega0;
  T omega_new = 1/(tau_t+tau_0);
  return omega_new;
}


///////////////////SHEAR IMPROVED SMAGORINSKY//////////////////////////
template<typename T, template<typename U> class Lattice>
ShearSmagorinskyBGKdynamics<T,Lattice>::ShearSmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_), preFactor(smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2)), dx(dx_), dt(dt_)
{ }

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyBGKdynamics<T,Lattice>::collide( Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), rho, pi, cell, statistics.getTime());
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), rho, pi, cell, statistics.getTime());
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyBGKdynamics<T,Lattice>::setOmega(T omega_)
{
//  this->omega = omega_;
  BGKdynamics<T,Lattice>::setOmega(omega_);
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell,
    int iT)
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), rho, pi, cell, iT);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyBGKdynamics<T,Lattice>::computeOmega(T omega0,
    T rho, T pi[util::TensorVal<Lattice<T> >::n], Cell<T,Lattice>& cell, int iT )
{
  OstreamManager clout(std::cout,"shearImprovedCollide");

  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.getExternal(avShearIsAt);
  *avShear = (*avShear*iT+PiNeqNorm)/(iT+1);
  //clout << "avShear"<< *avShear<<endl;
  T tau_0 = 1./omega0;
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(preFactor*PiNeqNorm_SISM))-tau_0);

  T omega_new = 1./(tau_t+tau_0);
  //clout << iT << std::endl;
  return omega_new;
}


///////////////////////// FORCED SHEAR SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
ShearSmagorinskyForcedBGKdynamics<T,Lattice>::ShearSmagorinskyForcedBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : ForcedBGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_), preFactor(smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2)), dx(dx_), dt(dt_)
{ }

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyForcedBGKdynamics<T,Lattice>::collide( Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), rho, pi, cell, statistics.getTime());

  T* force = cell.getExternal(this->forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}



template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyForcedBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  ForcedBGKdynamics<T,Lattice>::setOmega(omega_);
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyForcedBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell,
    int iT)
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), rho, pi, cell, iT);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyForcedBGKdynamics<T,Lattice>::computeOmega(T omega0,
    T rho, T pi[util::TensorVal<Lattice<T> >::n], Cell<T,Lattice>& cell, int iT )
{
  OstreamManager clout(std::cout,"shearImprovedCollide");

  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.getExternal(avShearIsAt);
  *avShear = (*avShear*iT+PiNeqNorm)/(iT+1);
  //clout << "avShear"<< *avShear<<endl;

  T tau_0 = 1./omega0;
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(preFactor*PiNeqNorm_SISM))-tau_0);

  T omega_new = 1./(tau_t+tau_0);
  //clout << iT << std::endl;
  return omega_new;
}



////////////////////// Class SmagorinskyBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
SmagorinskyBGKdynamics<T,Lattice>::SmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyBGKdynamics<T,Lattice>::staticCollide(Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyBGKdynamics<T,Lattice>::setOmega(T omega_)
{
//  this->omega = omega_;
  BGKdynamics<T,Lattice>::setOmega(omega_);
  preFactor = computePreFactor(omega_, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}



template<typename T, template<typename U> class Lattice>
T SmagorinskyBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor_, T rho,
    T pi[util::TensorVal<Lattice<T> >::n] )
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;

}




///////////////////////// FORCED SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
SmagorinskyForcedBGKdynamics<T,Lattice>::SmagorinskyForcedBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : ForcedBGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyForcedBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T* force = cell.getExternal(this->forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyForcedBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyForcedBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  ForcedBGKdynamics<T,Lattice>::setOmega(omega_);
  preFactor = computePreFactor(omega_, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyForcedBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyForcedBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}



template<typename T, template<typename U> class Lattice>
T SmagorinskyForcedBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor_, T rho,
    T pi[util::TensorVal<Lattice<T> >::n] )
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;

}

///////////////////////// FORCED Linear Velocity SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::SmagorinskyLinearVelocityForcedBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : ForcedBGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T* force = cell.getExternal(this->forceBeginsAt);
  int nDim = Lattice<T>::d;
  T forceSave[nDim];
  // adds a+Bu to force, where
  //   d=2: a1=v[0], a2=v[1], B11=v[2], B12=v[3], B21=v[4], B22=v[5]
  //   d=2: a1=v[0], a2=v[1], a3=v[2], B11=v[3], B12=v[4], B13=v[5], B21=v[6], B22=v[7], B23=v[8], B31=v[9], B32=v[10], B33=v[11]
  T* v = cell.getExternal(Lattice<T>::ExternalField::vBeginsAt);
  for (int iDim=0; iDim<nDim; ++iDim) {
    forceSave[iDim] = force[iDim];
    force[iDim] += v[iDim];
    for (int jDim=0; jDim<nDim; ++jDim) {
      force[iDim] += v[jDim + iDim*nDim + nDim]*u[jDim];
    }
  }
  for (int iVel=0; iVel<nDim; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
  // Writing back to froce fector
  for (int iVel=0; iVel<nDim; ++iVel) {
    force[iVel] = forceSave[iVel];
  }
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  ForcedBGKdynamics<T,Lattice>::setOmega(omega_);
  preFactor = computePreFactor(omega_, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}



template<typename T, template<typename U> class Lattice>
T SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor_, T rho,
    T pi[util::TensorVal<Lattice<T> >::n] )
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;

}

////////////////////// Class KrauseBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
KrauseBGKdynamics<T,Lattice>::KrauseBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_)
  : BGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void KrauseBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  T newOmega[Lattice<T>::q];
  this->_momenta.computeRhoU(cell, rho, u);
  computeOmega(this->getOmega(), cell, preFactor, rho, u, newOmega);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void KrauseBGKdynamics<T,Lattice>::staticCollide(Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  T newOmega[Lattice<T>::q];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  computeOmega(this->getOmega(), cell, preFactor, rho, uTemp, newOmega);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  if (cell.takesStatistics()) {
    statistics.incrementStats(rho, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void KrauseBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  BGKdynamics<T,Lattice>::setOmega(omega_);
  preFactor = computePreFactor(omega_, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T KrauseBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*3*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}



template<typename T, template<typename U> class Lattice>
void KrauseBGKdynamics<T,Lattice>::computeOmega(T omega0, Cell<T,Lattice>& cell, T preFactor_, T rho,
    T u[Lattice<T>::d], T newOmega[Lattice<T>::q])
{
  T uSqr = u[0]*u[0];
  for (int iDim=0; iDim<Lattice<T>::d; iDim++) {
    uSqr += u[iDim]*u[iDim];
  }
  /// Molecular realaxation time
  T tau_mol = 1./omega0;

  for (int iPop=0; iPop<Lattice<T>::q; iPop++) {
    T fNeq = std::fabs(cell[iPop] - lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr));
    /// Turbulent realaxation time
    T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*fNeq) - tau_mol);
    /// Effective realaxation time
    tau_eff = tau_mol + tau_turb;
    newOmega[iPop] = 1./tau_eff;
  }
}



}

#endif
