/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2011 Jonas Latt, Mathias J. Krause,
 *  Jonas Kratzke
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

#ifndef UNITS_HH
#define UNITS_HH

#include "units.h"

namespace olb {

template<typename T>
void writeLogFile(LBunits<T> const& converter,
                  std::string const& title)
{
  std::string fullName = singleton::directories().getLogOutDir() + "olbLog.dat";
  olb_ofstream ofile(fullName.c_str());
  ofile << title << "\n\n";
  ofile << "Velocity in lattice units: u="         << converter.getLatticeU()   << "\n";
  ofile << "Reynolds number:           Re="        << converter.getRe()         << "\n";
  ofile << "Lattice viscosity:         latticeNu=" << converter.getLatticeNu()  << "\n";
  ofile << "Lattice resolution:        N="         << converter.getResolution() << "\n";
  ofile << "Extent of the system:      lx="        << converter.getLx()         << "\n";
  ofile << "Extent of the system:      ly="        << converter.getLy()         << "\n";
  ofile << "Extent of the system:      lz="        << converter.getLz()         << "\n";
  ofile << "Grid spacing deltaX:       dx="        << converter.getDeltaX()     << "\n";
  ofile << "Time step deltaT:          dt="        << converter.getDeltaT()     << "\n";
  ofile << "Relaxation time:           tau="       << converter.getTau()        << "\n";
  ofile << "Relaxation frequency:      omega="     << converter.getOmega()      << "\n";
}

template<typename T>
void writeLogFile(LBconverter<T> const& converter,
                  std::string const& title)
{
  std::string fullName = singleton::directories().getLogOutDir() + title + ".dat";
  olb_ofstream ofile(fullName.c_str());
  ofile << "LBconverter information\n\n";
  ofile << "characteristical values\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "Dimension(d):                     dim="              << converter.getDim()              << "\n";
  ofile << "Characteristical length(m):       charL="            << converter.getCharL()            << "\n";
  ofile << "Characteristical speed(m/s):      charU="            << converter.getCharU()            << "\n";
  ofile << "Characteristical time(s):         charT="            << converter.getCharTime()         << "\n";
  ofile << "Density factor(kg/m^d):           charRho="          << converter.getCharRho()          << "\n";
  ofile << "Characterestical mass(kg):        charMass="         << converter.getCharMass()         << "\n";
  ofile << "Characterestical force(N):        charForce="        << converter.getCharForce()        << "\n";
  ofile << "Characterestical pressure(Pa):    charPressure="     << converter.getCharPressure()     << "\n";
  ofile << "Pressure level(Pa):               pressureLevel="    << converter.getPressureLevel()    << "\n";
  ofile << "Phys. kinematic viscosity(m^2/s): charNu="           << converter.getCharNu()           << "\n";
  ofile << "Phys. dynamic viscosity(N*s/m^2): dynVisco="         << converter.getDynamicViscosity() << "\n";
  ofile << "======================================================================\n\n";
  ofile << "lattice values\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "DeltaX:                           deltaX="           << converter.getDeltaX()           << "\n";
  ofile << "Lattice velocity:                 latticeU="         << converter.getLatticeU()         << "\n";
  ofile << "DeltaT:                           deltaT="           << converter.getDeltaT()           << "\n";
  ofile << "Reynolds number:                  Re="               << converter.getRe()               << "\n";
  ofile << "DimlessNu:                        dNu="              << converter.getDimlessNu()        << "\n";
  ofile << "Viscosity for computation:        latticeNu="        << converter.getLatticeNu()        << "\n";
  ofile << "Relaxation time:                  tau="              << converter.getTau()              << "\n";
  ofile << "Relaxation frequency:             omega="            << converter.getOmega()            << "\n";
  ofile << "======================================================================\n\n";
  ofile << "conversion factors\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "latticeL(m):                      latticeL="         << converter.getLatticeL()         << "\n";
  ofile << "Time step (s):                    physTime="         << converter.physTime()            << "\n";
  ofile << "Velocity factor(m/s):             physVelocity="     << converter.physVelocity()        << "\n";
  ofile << "FlowRate factor(m^d/s):           physFlowRate="     << converter.physFlowRate()        << "\n";
  ofile << "Mass factor(kg):                  physMass="         << converter.physMass()            << "\n";
  ofile << "Force factor(N):                  physForce="        << converter.physForce()           << "\n";
  ofile << "Force factor massless(N/kg):      physMasslessForce="<< converter.physMasslessForce()   << "\n";
  ofile << "Pressure factor(Pa):              physPressure="     << converter.physPressure()        << "\n";
  ofile << "latticePressure:                  latticeP="         << converter.latticePressure()     << "\n";

}

template<typename T>
void LBconverter<T>::print() const
{
  clout << "LBconverter information" << std::endl;
  clout << "characteristical values" << std::endl;
  clout << "Dimension(d):                     dim="              << getDim()              << std::endl;
  clout << "Characteristical length(m):       charL="            << getCharL()            << std::endl;
  clout << "Characteristical speed(m/s):      charU="            << getCharU()            << std::endl;
  clout << "Characteristical time(s):         charT="            << getCharTime()         << std::endl;
  clout << "Density factor(kg/m^d):           charRho="          << getCharRho()          << std::endl;
  clout << "Characterestical mass(kg):        charMass="         << getCharMass()         << std::endl;
  clout << "Characterestical force(N):        charForce="        << getCharForce()        << std::endl;
  clout << "Characterestical pressure(Pa):    charPressure="     << getCharPressure()     << std::endl;
  clout << "Pressure level(Pa):               pressureLevel="    << getPressureLevel()    << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s): charNu="           << getCharNu()           << std::endl;

  clout << "lattice values" << std::endl;
  clout << "DeltaX:                           deltaX="           << getDeltaX()           << std::endl;
  clout << "Lattice velocity:                 latticeU="         << getLatticeU()         << std::endl;
  clout << "DeltaT:                           deltaT="           << getDeltaT()           << std::endl;
  clout << "Reynolds number:                  Re="               << getRe()               << std::endl;
  clout << "DimlessNu:                        dNu="              << getDimlessNu()        << std::endl;
  clout << "Viscosity for computation:        latticeNu="        << getLatticeNu()        << std::endl;
  clout << "Relaxation time:                  tau="              << getTau()              << std::endl;
  clout << "Relaxation frequency:             omega="            << getOmega()            << std::endl;

  clout << "conversion factors" << std::endl;
  clout << "latticeL(m):                      latticeL="         << getLatticeL()         << std::endl;
  clout << "Time step (s):                    physTime="         << physTime()            << std::endl;
  clout << "Velocity factor(m/s):             physVelocity="     << physVelocity()        << std::endl;
  clout << "FlowRate factor(m^d/s):           physFlowRate="     << physFlowRate()        << std::endl;
  clout << "Mass factor(kg):                  physMass="         << physMass()            << std::endl;
  clout << "Force factor(N):                  physForce="        << physForce()           << std::endl;
  clout << "Force factor massless(N/kg):      physMasslessForce="<< physMasslessForce()   << std::endl;
  clout << "Pressure factor(Pa):              physPressure="     << physPressure(4)       << std::endl;
  clout << "latticePressure:                  latticeP="         << latticePressure()     << std::endl;

}

template<typename T>
LBconverter<T>* createLBconverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createLBconverter");

  int dim = 0;
  T latticeL = T(0);
  T deltaX = T(0);
  int N = 0;
  T latticeU = T(0);
  T charNu = T(0);
  T Re = T(0);
  T charL = T(1);
  T charU = T(1);
  T charRho = T(1);
  T pressureLevel = T(0);
  bool verbose = false;
  params.setWarningsOn(false);
  // params[parameter].read(value) sets the value or returns false if the parameter can not be found

  if (!params["Application"]["dim"].read<int>(dim, verbose)) {
    clout << "Error: Cannot read parameter from Xml-file: dim" << std::endl;
    exit (1);
  }
  if (!params["Application"]["DiscParam"]["latticeU"].read(latticeU, verbose)) {
    clout << "Error: Cannot read parameter from Xml-file: latticeU"
          << std::endl;
    exit (1);
  }
  if (!params["Application"]["PhysParam"]["charL"].read(charL, verbose)) {
    clout << "Parameter charL not found in Xml-file. Set default: charL = 1."
          << std::endl;
  }
  if (!params["Application"]["DiscParam"]["latticeL"].read(latticeL, verbose)) {
    if (!params["Application"]["DiscParam"]["deltaX"].read(deltaX, verbose)) {
      if (!params["Application"]["DiscParam"]["resolution"].read<int>(N, verbose)) {
        clout << "Error: Cannot read any of the parameters from Xml-file: "
              << "latticeL, deltaX, resolution"
              << std::endl;
        exit (1);
      }
      deltaX = (T) 1 / (T) N;
    }
    latticeL = deltaX * charL;
  }
  if (!params["Application"]["PhysParam"]["charU"].read(charU, verbose)) {
    clout << "Parameter charU not found in Xml-file. Set default: charU = 1."
          << std::endl;
  }
  if (!params["Application"]["PhysParam"]["charRho"].read(charRho, verbose)) {
    clout << "Parameter charRho not found in Xml-file. Set default: charRho = 1."
          << std::endl;
  }
  if (!params["Application"]["PhysParam"]["charPressure"].read(pressureLevel, verbose)) {
    clout << "Parameter charPressure not found in Xml-file. Set default: charPressure = 0."
          << std::endl;
  }
  if (!params["Application"]["PhysParam"]["charNu"].read(charNu, verbose)) {
    if (!params["Application"]["PhysParam"]["Re"].read(Re, verbose)) {
      clout << "Error: Cannot read neither charNu nor Re from Xml-file."
            << std::endl;
      exit (1);
    }
    charNu = charL*charU / Re;
  }
  params.setWarningsOn(true);

  return new LBconverter<T>(dim, latticeL, latticeU, charNu, charL, charU,
                            charRho, pressureLevel);
}


////// ctor
template<typename T>
RTLBConstconverter<T>::RTLBConstconverter
(T latticeL, T physAbsorption, T physScattering, T physRefrMedia, T physRefrAmbient)
  : LBconverter<T>(3, latticeL, latticeL, 1./6. ), _physAbsorption(physAbsorption),
    _physScattering(physScattering), _physRefrMedia(physRefrMedia),
    _physRefrAmbient(physRefrAmbient), clout(std::cout,"RTLBconverter")
{
  _sinkTerm = 3 *_physAbsorption *(_physAbsorption+_physScattering) / 8. * latticeL * latticeL;
  computeZeta();
}

// public methods
template<typename T>
T RTLBConstconverter<T>::getAbsorption() const
{
  return _physAbsorption;
}

template<typename T>
T RTLBConstconverter<T>::getScattering() const
{
  return _physScattering;
}

template<typename T>
T RTLBConstconverter<T>::getSingleScatAlbedo() const
{
  return _physScattering / (_physAbsorption +_physScattering);
}

template<typename T>
T RTLBConstconverter<T>::getExtinctionCoeff() const
{
  return (_physAbsorption +_physScattering) * this->getLatticeL();
}

template<typename T>
T RTLBConstconverter<T>::getSinkTerm() const
{
  return _sinkTerm;
}
template<typename T>
T RTLBConstconverter<T>::getZeta() const
{
  return _latticeZeta;
}
template<typename T>
T RTLBConstconverter<T>::getRefractionAmbient() const
{
  return _physRefrAmbient;
}
template<typename T>
T RTLBConstconverter<T>::getRefractionMedia() const
{
  return _physRefrMedia;
}


template<typename T>
void RTLBConstconverter<T>::setZeta(  T zeta )
{
  _latticeZeta = zeta;
}


template<typename T>
void RTLBConstconverter<T>::print() const
{
  clout << "RTLBConstconverter information" << std::endl;
  clout << "pysical values" << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s): charNu="           << this->getCharNu()           << std::endl;

  clout << "lattice values" << std::endl;
  clout << "DeltaX:                           deltaX= latticeL/charL=" << this->getDeltaX()     << std::endl;
  clout << "DeltaT:                           deltaT="           << this->getDeltaT()           << std::endl;
  clout << "DimlessNu:                        dNu="              << this->getDimlessNu()        << std::endl;
  clout << "Viscosity for computation:        latticeNu="        << this->getLatticeNu()        << std::endl;
  clout << "Relaxation time:                  tau="              << this->getTau()              << std::endl;
  clout << "Relaxation frequency:             omega="            << this->getOmega()            << std::endl;
  clout << "Absorption:                       absorption 1/m="   << getAbsorption()             << std::endl;
  clout << "Scattering:                       scattering 1/m="   << getScattering()             << std::endl;
  clout << "Extinction:                       extinction= "      << _physAbsorption +_physScattering  << std::endl;
  clout << "ScatAlbedo:                       scatAlbedo= "      << this->getSingleScatAlbedo() << std::endl;

  clout << "Mesoscopic sink term              sink="             << getSinkTerm()               << std::endl;
  clout << "Zeta partial BB:                  zeta partial BB="  << getZeta()                   << std::endl;


  clout << "conversion factors" << std::endl;
  clout << "latticeL(m):                      latticeL="         << this->getLatticeL()         << std::endl;
  clout << "Time step (s):                    physTime="         << this->physTime()            << std::endl;

}

// privat methods
template<typename T>
T RTLBConstconverter<T>::theta_ (T theta)
{
  T nrel = _physRefrMedia / _physRefrAmbient;
  T theta_ = asin( nrel * sin(theta));
  return theta_;
}

template<typename T>
T RTLBConstconverter<T>::rf (T theta)
{
  T nrel = _physRefrMedia / _physRefrAmbient;
  T rf_1 = 0.5 * pow((nrel * cos(theta_(theta)) - cos(theta)) /
                     (nrel * cos(theta_(theta)) + cos(theta)), 2.);
  T rf_2 = 0.5 * pow((nrel * cos(theta) - cos(theta_(theta))) /
                     (nrel * cos(theta) + cos(theta_(theta))), 2.);
  T rf = rf_1 + rf_2;
  return rf;
}

template<typename T>
T RTLBConstconverter<T>::r_phi_diff (T theta)
{
  T r_phi_diff = 2. * sin(theta) * cos(theta) * rf(theta);
  return r_phi_diff;
}

template<typename T>
T RTLBConstconverter<T>::r_j_diff (T theta)
{
  T r_j_diff = 3. * sin(theta) * pow(cos(theta),2.) * rf(theta);
  return r_j_diff;
}

template<typename T>
void RTLBConstconverter<T>::computeZeta ()
{
  T n = 10000.0;
  T h = (M_PI / 2.) /n;
  T r_phi = 0.0;
  T r_j = 0.0;
  for (int i = 0; i < n; i++) {
    r_phi += h*(r_phi_diff(0.5*h + h*i));
    r_j   += h*(r_j_diff  (0.5*h + h*i));
  }
  T r_eff = (r_phi + r_j) / (2 - r_phi + r_j);
  T c_r  = (1 + r_eff) / (1 - r_eff);
  T epsilon  = 1.5 * (_physAbsorption + _physScattering) * this->getLatticeL();
  _latticeZeta = 2. / (1. + (epsilon / (2. * c_r)));
}

// write logfile
template<typename T>
void writeLogFile(RTLBConstconverter<T> const& converter, std::string const& title)
{
  std::string fullName = singleton::directories().getLogOutDir() + title + ".dat";
  olb_ofstream ofile(fullName.c_str());
  ofile << "LBconverter information\n\n";
  ofile << "physical values\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "Dimension(d):                     dim="              << converter.getDim()              << "\n";
  ofile << "Phys. kinematic viscosity(m^2/s): charNu="           << converter.getCharNu()           << "\n";
  ofile << "Phys. dynamic viscosity(N*s/m^2): dynVisco="         << converter.getDynamicViscosity() << "\n";
  ofile << "Absorption(m^-1):                 absorption="       << converter.getAbsorption()       << "\n";
  ofile << "Scattering(m^-1):                 scattering="       << converter.getScattering()       << "\n";

  ofile << "======================================================================\n\n";
  ofile << "lattice values\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "DeltaX:                           deltaX="           << converter.getDeltaX()           << "\n";
  ofile << "Lattice velocity:                 latticeU="         << converter.getLatticeU()         << "\n";
  ofile << "DeltaT:                           deltaT="           << converter.getDeltaT()           << "\n";
  ofile << "Reynolds number:                  Re="               << converter.getRe()               << "\n";
  ofile << "DimlessNu:                        dNu="              << converter.getDimlessNu()        << "\n";
  ofile << "Viscosity for computation:        latticeNu="        << converter.getLatticeNu()        << "\n";
  ofile << "Relaxation time:                  tau="              << converter.getTau()              << "\n";
  ofile << "Relaxation frequency:             omega="            << converter.getOmega()            << "\n";
  ofile << "Sink term:                        eta="              << converter.getSinkTerm()         << "\n";
  ofile << "Knutsen number:                   kappa= "           << 1/(converter.getExtinctionCoeff()/converter.getLatticeL()) << "\n";
  //ofile << "Extinction:                       extinction= "      << _physAbsorption +_physScattering()  << "\n";
  ofile << "ScatAlbedo:                       scatAlbedo= "      << converter.getSingleScatAlbedo() << "\n";
  ofile << "======================================================================\n\n";
  ofile << "conversion factors\n";
  ofile << "----------------------------------------------------------------------\n";
  ofile << "latticeL(m):                      latticeL="         << converter.getLatticeL()         << "\n";
  ofile << "Time step (s):                    physTime="         << converter.physTime()            << "\n";

}


}  // namespace olb

#endif
