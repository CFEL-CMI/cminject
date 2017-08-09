/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

#ifndef ADVECTION_DIFFUSION_UNITS_H
#define ADVECTION_DIFFUSION_UNITS_H

#include <string>
#include <fstream>
#include "core/singleton.h"

namespace olb {

/// A useful class for the conversion between dimensionless and lattice units.
template<typename T, template<typename NSU> class NSLattice, template<typename ADU> class ADLattice>
class AdvectionDiffusionUnitLB {
public:
  /// Constructor
  /** \param Ra               Raylegh number
   *  \param Pr               Prandtl number
   *  \param Tcold            minimum temperature
   *  \param Thot             maximum temperature
   *  \param deltaTemperature Reynolds number
   *  \param deltaT           time discretization number
   *  \param N                resolution (a lattice of size 1 has N_+1 cells)
   *  \param lx               x-length in dimensionless units (e.g. 1)
   *  \param ly               y-length in dimensionless units (e.g. 1)
   *  \param lz               z-length in dimensionless units (e.g. 1)
   */
  AdvectionDiffusionUnitLB(T Ra, T Pr, T Tcold, T Thot, T N, T deltaT,
                           T lx, T ly, T lz=T() )
    : _Ra(Ra), _Pr(Pr), _Tcold(Tcold), _Thot(Thot), _NN(N), _deltaT(deltaT),
      _lx(lx), _ly(ly), _lz(lz)
  {
  }
  /// Rayleigh number
  T getRa() const
  {
    return _Ra;
  }
  /// Prandlt number
  T getPr() const
  {
    return _Pr;
  }
  /// minimum temperature
  T getTcold() const
  {
    return _Tcold;
  }
  /// maximum temperature
  T getThot() const
  {
    return _Thot;
  }
  /// delta temperature number
  T getDeltaTemperature() const
  {
    return _Thot - _Tcold;
  }
  /// mid-temperature
  T getT0() const
  {
    return (_Thot + _Tcold)/(T)2;
  }
  /// resolution (a lattice of size 1 has getN()+1 cells)
  T getN() const
  {
    return _NN;
  }
  /// x-length in dimensionless units
  T getLx() const
  {
    return _lx;
  }
  /// y-length in dimensionless units
  T getLy() const
  {
    return _ly;
  }
  /// z-length in dimensionless units
  T getLz() const
  {
    return _lz;
  }
  /// lattice spacing in dimensionless units
  T getDeltaX() const
  {
    return (T)1 / _NN;
  }
  /// time step in dimensionless units
  T getDeltaT() const
  {
    return _deltaT;
  }
  /// conversion from dimensionless to lattice units for space coordinate
  int nCell(T l) const
  {
    return (int)(l/getDeltaX()+(T)0.5);
  }
  /// conversion from dimensionless to lattice units for time coordinate
  int nStep(T t) const
  {
    return (int)(t/getDeltaT()+(T)0.5);
  }
  /// number of lattice cells in x-direction
  int getNx() const
  {
    return nCell(_lx) + 1;
  }
  /// number of lattice cells in y-direction
  int getNy() const
  {
    return nCell(_ly) + 1;
  }
  /// number of lattice cells in z-direction
  int getNz() const
  {
    return nCell(_lz) + 1;
  }
  /// velocity in lattice units (proportional to Mach number)
  T getU() const
  {
    return getDeltaT() / getDeltaX()  ;
  }
  /// viscosity in lattice units
  T getNu() const
  {
    return sqrt(getPr()/getRa()) * getDeltaT() / (getDeltaX()*getDeltaX());
  }
  /// thermal conductivity in lattice units
  T getKappa() const
  {
    return sqrt((T)1/(getPr()*getRa()))*getDeltaT()/(getDeltaX()*getDeltaX());
  }
  /// viscosity in lattice units
  T getGravity() const
  {
    return getDeltaT() * getDeltaT() / getDeltaX();
  }
  /// relaxation time
  T getTauNS() const
  {
    return NSLattice<T>::invCs2*getNu()+(T)0.5;
  }
  /// relaxation frequency
  T getOmegaNS() const
  {
    return (T)1 / getTauNS();
  }
  /// relaxation time
  T getTauT() const
  {
    return ADLattice<T>::invCs2*getKappa() + (T)0.5;
  }
  /// relaxation frequency
  T getOmegaT() const
  {
    return (T)1 / getTauT();
  }
private:
  T _Ra, _Pr, _Tcold, _Thot, _NN, _deltaT, _lx, _ly, _lz;
};

template<typename T, template<typename NSU> class NSLattice, template<typename ADU> class ADLattice>
void writeLogFile(AdvectionDiffusionUnitLB<T,NSLattice,ADLattice> const& converter,
                  std::string const& title)
{
  std::string fullName = singleton::directories().getLogOutDir() + "olbLog.dat";
  std::ofstream ofile(fullName.c_str());
  ofile << title << "\n\n";
  ofile << "Velocity in lattice units: u="  << converter.getU() << "\n";
  ofile << "Raynleigh number:          Ra=" << converter.getRa() << "\n";
  ofile << "Prandlt number:            Pr=" << converter.getPr() << "\n";
  ofile << "Kinematic viscosity:       Nu=" << converter.getNu() << "\n";
  ofile << "AdvectionDiffusion conductivity:      Kappa=" << converter.getKappa() << "\n";
  ofile << "Lattice resolution:        N="  << converter.getN() << "\n";
  ofile << "Extent of the system:      lx=" << converter.getLx() << "\n";
  ofile << "Extent of the system:      ly=" << converter.getLy() << "\n";
  ofile << "Extent of the system:      lz=" << converter.getLz() << "\n";
  ofile << "Grid spacing deltaX:       dx=" << converter.getDeltaX() << "\n";
  ofile << "Time step deltaT:          dt=" << converter.getDeltaT() << "\n";
}

}  // namespace olb

#endif
