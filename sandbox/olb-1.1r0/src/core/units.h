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

/** \file
 * Unit handling -- header file.
 */

#ifndef UNITS_H
#define UNITS_H

//#include "communication/mpiManager.h"
#include "io/parallelIO.h"
#include "io/xmlReader.h"
#include "io/ostreamManager.h"
#include <string>
#include <fstream>
#include "singleton.h"
#include "olbDebug.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/// All OpenLB code is contained in this namespace.
namespace olb {

/// Conversion between dimensionless and lattice units with on-lattice boundaries
template<typename T>
class LBunits {
public:
  /** Constructor:
   *  \param latticeU velocity in lattice units (proportional to Mach number)
   *  \param Re       Reynolds number
   *  \param N        resolution (a lattice of size 1 has N_+1 cells)
   *  \param lx       x-length in dimensionless units (e.g. 1)
   *  \param ly       y-length in dimensionless units (e.g. 1)
   *  \param lz       z-length in dimensionless units (e.g. 1)
   */
  LBunits(T latticeU, T Re, int resolution, T lx, T ly, T lz=T() )
    : _latticeU(latticeU), _Re(Re), _resolution(resolution), _lx(lx), _ly(ly), _lz(lz)
  {
  }
  /// velocity in lattice units (proportional to Mach number)
  T getLatticeU() const
  {
    return _latticeU;
  }
  /// Reynolds number
  T getRe() const
  {
    return _Re;
  }
  /// resolution
  int getResolution() const
  {
    return _resolution;
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
    return (T)1/(T)getResolution();
  }
  /// time step in dimensionless units
  T getDeltaT() const
  {
    return getDeltaX()*getLatticeU();
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
  int getNx(bool offLattice=false) const
  {
    return nCell(_lx)+1+(int)offLattice;
  }
  /// number of lattice cells in y-direction
  int getNy(bool offLattice=false) const
  {
    return nCell(_ly)+1+(int)offLattice;
  }
  /// number of lattice cells in z-direction
  int getNz(bool offLattice=false) const
  {
    return nCell(_lz)+1+(int)offLattice;
  }
  /// viscosity in lattice units
  T getLatticeNu() const
  {
    return getLatticeU()*getResolution()/_Re;
  }
  /// relaxation time
  T getTau() const
  {
    return (T)3*getLatticeNu()+(T)0.5;
  }
  /// relaxation frequency
  T getOmega() const
  {
    return (T)1 / getTau();
  }
private:
  T _latticeU;
  T _Re;
  int _resolution;
  T _lx;
  T _ly;
  T _lz;
};

template<typename T>
void writeLogFile(LBunits<T> const& converter, std::string const& title);

/// Conversion between lattice units and physical units, physXYZ(lattice quantity XYZ) returns the physical size of a lattice quantity XYZ
template<typename T>
class LBconverter {
public:
  /** Constructor:
   *  \param dim            dimension of the domain (2D or 3D)
   *  \param latticeL       length of a lattice cell in meter (proportional to Knudsen number)
   *  \param latticeU       velocity in dimensionless lattice units (proportional to Mach number)
   *  \param charNu         kinematic viscosity in m^2/s
   *  \param charL          characteristical length in meter
   *  \param charU          characteristical speed in m/s
   *  \param charRho        density factor in kg/m^d (latticeRho can be multplied by this factor
   *                        to get the local physical density)
   *  \param pressureLevel additive pressure constant in Pa (added to the relative pressure
   *                       result of the computation to get the absolute value)
   */
  LBconverter(int dim, T latticeL, T latticeU, T charNu, T charL = 1, T charU = 1,
              T charRho = 1, T pressureLevel = 0 )
    : clout(std::cout,"LBconverter"), _dim(dim), _latticeL(latticeL), _latticeU(latticeU),
      _charNu(charNu), _charL(charL), _charU(charU), _charRho(charRho),
      _pressureLevel(pressureLevel)
  {
    singleton::checkValue(_charNu);
  }
  // virtual destructor
  virtual ~LBconverter() {};
  /// dimension of the domain (2D or 3D)
  int getDim() const
  {
    return _dim;
  }
  /// length of a lattice cell in meter
  T getLatticeL() const
  {
    singleton::checkValue(_charNu);
    return _latticeL;
  }
  /// characteristical length in meter
  T getCharL(T dL = 1) const
  {
    return _charL * dL;
  }
  /// characteristical speed in m/s
  T getCharU(T dU = 1) const
  {
    return _charU * dU;
  }
  /// characteristical time in s
  T getCharTime(T dT = 1) const
  {
    return _charL/_charU * dT;
  }
  /// kinematic viscosity in m^2/s
  T getCharNu() const
  {
    return _charNu;
  }
  /// dynamic viscosity in N*s/m^2
  T getDynamicViscosity() const
  {
    return _charNu * _charRho;
  }
  /// density factor in kg/m^d
  T getCharRho(T dRho = 1) const
  {
    return _charRho * dRho;
  }
  /// characteristical mass in kg
  T getCharMass(T dM = 1) const
  {
    return _charRho*pow(_charL,_dim) * dM;
  }
  /// characteristical force in Newton = kg*m/s^2
  T getCharForce(T dF = 1) const
  {
    return getCharMass()*_charL / (getCharTime()*getCharTime()) * dF;
  }
  /// characteristical pressure in Pascal = N/m^(d-1) = rho*m^2/t^2
  T getCharPressure(T dP = 1) const
  {
    return getCharForce() / pow(_charL, _dim-1 ) * dP + _pressureLevel;
  }
  /// characteristical Pressure in Pa
  T getPressureLevel() const
  {
    return _pressureLevel;
  }

  /// Reynolds number
  T getRe() const
  {
    return _charL * _charU / _charNu;
  }
  /// dimensionless kinematic viscosity
  T getDimlessNu() const
  {
    return 1 / getRe();
  }
  /// discretization parameter for grid-spacing (proportional to Knudsen number)
  T getDeltaX() const
  {
    return _latticeL / _charL;
  }
  /// discretization parameter for velocity (proportional to Mach number)
  T getLatticeU() const
  {
    return _latticeU;
  }
  /// discretization parameter for time
  T getDeltaT() const
  {
    return _latticeU * getDeltaX();
  }
  /// lattice kinematic viscosity used for computation
  T getLatticeNu() const
  {
    return getDeltaT()/ (getDeltaX() * getDeltaX() * getRe());
  }
  /// relaxation time
  T getTau() const
  {
    return (T)3*getLatticeNu()+(T)0.5;
  }
  /// relaxation frequency
  T getOmega() const
  {
    return (T)1 / getTau();
  }

  /// physical length of a number of cells
  T physLength(T latticeLength = 1) const
  {
    return _charL * getDeltaX() * latticeLength;
  }
  /// length of a lattice time period in seconds
  /// default: get conversion factor -> lattice to physical time
  T physTime(T latticeTime = 1) const
  {
    return _charL/_charU * getDeltaT() * latticeTime;
  }
  /// convert lattice velocity to physical velocity in m/s
  /// default: get conversion factor -> lattice to physical velocity
  T physVelocity(T latticeVelocity = 1) const
  {
    return _charU / _latticeU * latticeVelocity;
  }
  /// convert lattice flow rate to physical flow rate
  /// default: get conversion factor -> lattice to physical flow rate
  T physFlowRate(T latticeFlowRate = 1) const
  {
    return latticeFlowRate*pow(physLength(),_dim) / physTime();
  }
  /// convert lattice to physical density
  /// default: get conversion factor -> lattice to physical density
  T physRho(T latticeRho = 1) const
  {
    return _charRho*latticeRho;
  }
  /// convert lattice density to physical mass in kg
  /// default: get conversion factor -> lattice to physical mass
  T physMass(T latticeRho = 1) const
  {
    return physRho(latticeRho)*pow(physLength(),_dim);
  }
  /// convert lattice to physical force in Newton
  /// default: get conversion factor -> lattice to physical force
  T physForce(T latticeForce = 1) const
  {
    return physMass() * physLength() / (physTime() * physTime()) * latticeForce;
  }
  /// convert lattice to physical massless force in Newton/kg
  /// default: get conversion factor -> lattice to physical massless force
  T physMasslessForce(T latticeForce = 1) const
  {
    return physForce(latticeForce) / physMass();
  }
  /// convert: lattice to physical pressure in Pa
  /// physicalPressure = (rho-1)/3)*pressureFactor
  T physPressure(T latticePressure = 1) const
  {
    return latticePressure*physForce() / (pow(physLength(),_dim-1)) + _pressureLevel;
  }
  /// convert: lattice rho to physical pressure in Pa
  /// physicalPressure = (rho-1)/3)*pressureFactor
  T physPressureFromRho(T rho) const
  {
    return ((rho - 1) / 3.)*physForce() / (pow(physLength(),_dim-1)) + _pressureLevel;
  }

  /// convert: physical length to lattice length
  T latticeLength(T physicalLength = 1) const
  {
    return physicalLength / physLength();
  }
  /// convert: physical velocity to lattice velocity
  T latticeVelocity(T physicalVelocity = 1) const
  {
    return physicalVelocity / physVelocity();
  }
  /// returns number of lattice cells within a length l
  int numCells(T physicalLength = -1) const
  {
    singleton::checkValue(_charNu);
    if (physicalLength == -1) {
      physicalLength = _charL;
    }
    return (int)(physicalLength / physLength()+T(0.5) );
  }
  /// returns number of lattice nodes of a physical length l
  int numNodes(T physicalLength = -1) const
  {
    if (physicalLength == -1) {
      physicalLength = _charL;
    }
    return (int)(physicalLength / physLength()+(1.5));
  }
  /// returns number of lattice time steps within a period physicalT
  int numTimeSteps(T physicalTime) const
  {
    return (int)(physicalTime / physTime()+T(0.5));
  }
  /// convert physical to lattice pressure
  /// default: get conversion factor -> physical to lattice pressure
  T latticePressure(T physicalPressure = 1) const
  {
    return (physicalPressure - _pressureLevel) * (pow(physLength(),_dim-1))/physForce();
  }
  /// convert: physical pressure in Pa to lattice density
  /// latticeRho = physical pressure / pressureFactor * 3 -1
  T rhoFromPhysicalPressure(T physicalPressure = 0) const
  {
    return (physicalPressure - _pressureLevel) * (pow(physLength(),_dim-1))/physForce() *T(3) + T(1);
  }

  /// convert physical to lattice force
  /// default: get conversion factor -> physical to lattice force
  T latticeForce(T physicalForce = 1) const
  {
    return physicalForce / physForce();
  }

  /// converts a physical permeability K to a lattice-dependent porosity d
  /// (a velocity scaling factor depending on Maxwellian distribution function),
  /// needs PorousBGKdynamics
  T latticePorosity(T K) const
  {
    OLB_PRECONDITION(K >= pow( physLength(), getDim() - 1 ) * getLatticeNu() * getTau());
    return 1 - pow( physLength(), getDim() - 1 ) * getLatticeNu() * getTau() / K;
  }

  /// converts a lattice-dependent porosity d (a velocity scaling factor
  /// depending on Maxwellian distribution function) to a physical permeability K,
  /// needs PorousBGKdynamics
  T physPermeability(T d) const
  {
    return pow( physLength(), getDim() - 1 ) * getLatticeNu() * getTau() / ( 1 - d ) ;
  }

  T gridTerm() const
  {
    return pow( physLength(), getDim() - 1 ) * getLatticeNu() * getTau();
  }

  /// print converter information
  virtual void print() const;

private:
  mutable OstreamManager clout;
  int _dim;
  T _latticeL;
  T _latticeU;
  T _charNu;
  T _charL;
  T _charU;
  T _charRho;
  T _pressureLevel;
};

template<typename T>
void writeLogFile(LBconverter<T> const& converter, std::string const& title);

template<typename T>
LBconverter<T>* createLBconverter(XMLreader const& params);


/** Unit converter for reaction-diffusion equation, radiative transport.
 * Results presented at DSFD 2015 Edinburgh, D3Q7 only.
 *
 * 1. Macroscopic target equation reads
 * \f[ \Delta \Phi = \frac{\sigma_a}{D} \Phi  - Sink \f]
 * for diffusion coefficient \f$ D = \frac{1}{3(\sigma_a+\sigma_s)} \f$.
 *
 * 2. Correspondig kinetic equation (implementation at advectionDiffusionLbHelpers3D.hh)
 * \f[ \tilde f_i = f_i -\omega [ f_i -f_i^{eq} ] - \eta(\sigma_a,\sigma_s) f_i \f]
 *
 * 3. Relation to viscosity \f$ \nu \f$ is (indenpent of absorption and scattering) given by
 * \f[ \nu = \frac{1}{6} \f]
 *
 * 4. Relation to  mesoscopic sink term \f$ \eta \f$ is given by
 * \f[ \eta(\sigma_a,\sigma_s) =  \frac{3 \sigma_a(\sigma_a + \sigma_s)}{ 8N^2} \f]
 *
 * \param _physAbsorption    (RTE) absorption coefficient, \f$ \sigma_a \f$
 * \param _physScattering    (RTE) scattering coefficient, \f$ \sigma_s \f$
 * \param _sinkTerm          (mesoscopic) sink term
 * \param _physRefrAmbient    refractive index of ambient medium
 * \param _physRefrMedia      refractive index of scattering medium
 * \param _latticeZeta        mesoscopic boundary value
 *
** Expansion:
 * Documentation: Biomedical Optics, Lihong V. Wang Hsin-I Wu; Hiorth, Lad, Evje and Skjæveland 2008
 */
template<typename T>
class RTLBConstconverter final : public LBconverter<T> {
public:
  RTLBConstconverter(T latticeL, T physAbsorption, T physScattering, T physRefrMedia=1, T physRefrAmbient=1);
  T getAbsorption() const;
  T getScattering() const;
  T getSinkTerm() const;
  T getSingleScatAlbedo() const;
  T getExtinctionCoeff() const;
  T getZeta() const;
  T getRefractionAmbient() const;
  T getRefractionMedia() const;
  virtual void print() const override;
  void setZeta ( T zeta );

private:
  T _physAbsorption;
  T _physScattering;
  T _sinkTerm;
  T _physRefrMedia;
  T _physRefrAmbient;
  T _latticeZeta;
  T theta_ (T theta);
  T rf (T theta);
  T r_phi_diff (T theta);
  T r_j_diff (T theta);
  void computeZeta ();
  mutable OstreamManager clout;
};

template<typename T>
void writeLogFile(RTLBConstconverter<T> const& converter, std::string const& title);


}  // namespace olb

#endif
