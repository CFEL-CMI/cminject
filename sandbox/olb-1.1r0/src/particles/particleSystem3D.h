/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Thomas Henn
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

#ifndef PARTICLE_SYSTEM_3D_H
#define PARTICLE_SYSTEM_3D_H

#include <deque>
#include "contactDetection/contactDetection.h"
#include "geometry/superGeometry3D.h"
#include "core/units.h"
#include "core/blockLatticeStructure3D.h"
#include "forces/force3D.h"
#include "functors/analyticalF.h"
#include "boundaries/boundary3D.h"
#include "functors/frameChangeF3D.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticeInterpPhysVelocity3D;

template<typename T, template<typename U> class PARTICLETYPE>
class Force3D;
template<typename T, template<typename U> class PARTICLETYPE>
class Boundary3D;
template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;
template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSystem3D;
template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSysVtuWriter;
template<typename T>
class SuperParticleSysVtuWriterMag;

template<typename T, template<typename U> class PARTICLETYPE>
class SimulateParticles {
public:
  SimulateParticles(ParticleSystem3D<T, PARTICLETYPE>* ps)
    : _pSys(ps)
  {
  }
  inline void simulate(T dT)
  {
    _pSys->computeForce();
    _pSys->explicitEuler(dT);
    //_pSys->rungeKutta4(dT);
  }
private:
  ParticleSystem3D<T, PARTICLETYPE>* _pSys;

};


template<typename T>
class SimulateParticles<T, RotatingParticle3D> {
public:
  SimulateParticles(ParticleSystem3D<T, RotatingParticle3D>* ps)
    : _pSys(ps)
  {
  }
  inline void simulate(T dT)
  {
    _pSys->computeForce();
    _pSys->explicitEuler(dT);
    _pSys->integrateTorque(dT);
  }
private:
  ParticleSystem3D<T, RotatingParticle3D>* _pSys;
};

template<typename T>
class SimulateParticles<T, MagneticParticle3D> {
public:
  SimulateParticles(ParticleSystem3D<T, MagneticParticle3D>* ps)
    : _pSys(ps)
  {
  }
  inline void simulate(T dT)
  {
    _pSys->resetMag();
    _pSys->computeForce();
    _pSys->explicitEuler(dT);
    // calculates changes in orientation (dipole moment direction)
    // due to torque moments induced in
    // interpMagForceForMagP3D & magneticForceForMagP3D
    _pSys->integrateTorqueMag(dT);

  }
private:
  ParticleSystem3D<T, MagneticParticle3D>* _pSys;

};


template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D {
public:

  /// Default constructor for ParticleSystem
  ParticleSystem3D() = default;
  /// Constructor for ParticleSystem
  ParticleSystem3D(SuperGeometry3D<T>&, LBconverter<T>& conv);
  /// Copy constructor for ParticleSystem
  ParticleSystem3D(const ParticleSystem3D<T, PARTICLETYPE>& pS);
  /// Move constructor for ParticleSystem
  ParticleSystem3D(ParticleSystem3D<T, PARTICLETYPE> && pS);
  /// Destructor for ParticleSystem
  virtual ~ParticleSystem3D() {}

  virtual void simulate(T deltatime);

  /// Get number of particles in ParticleSystem
  int size();
  /// Get number of particles including shadow particles in ParticleSystem
  int sizeInclShadow() const;
  /// Get number of active particles in ParticleSystem
  int numOfActiveParticles();
  /// Get number of linked forces in ParticleSystem
  int numOfForces();
  /// Get number of particles in vincinity of material number mat
  int countMaterial(int mat = 1);

  /// Add a particle to ParticleSystem
  void addParticle(PARTICLETYPE<T>& p);
  /// Removes all particles from system
  void clearParticles();
  /// Add a force to ParticleSystem
  void addForce(std::shared_ptr<Force3D<T, PARTICLETYPE> > pF);
  /// Add a boundary to ParticleSystem
  void addBoundary(std::shared_ptr<Boundary3D<T, PARTICLETYPE> > pF);

  /// Get reference to a particle in the ParticleSystem
  PARTICLETYPE<T>& operator[](const int i);
  const PARTICLETYPE<T>& operator[](const int i) const;

  /// Set velocity of all particles to fluid velocity
  template<template<typename V> class DESCRIPTOR>
  void setVelToFluidVel(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> &);

  /// Set particle velocity to analytical velocity (e.g. as inital condition
  void setVelToAnalyticalVel(AnalyticalConst3D<T,T>&);

  /// Set global coordinates and extends of Particlesystem (SI units)
  void setPosExt(std::vector<T> physPos, std::vector<T> physExtend);

  /// Get global coordinates and extends of Particlesystem (SI units)
  const std::vector<T>& getPhysPos();
  const std::vector<T>& getPhysExtend();

  /// Save particle positions to file
  void saveToFile(std::string name);

  /// Compute all forces on particles
  void computeForce();
  /// Compute boundary contact
  void computeBoundary();

  /// Set boundary detection algorithm (for future features)
  void setContactDetection(ContactDetection<T, PARTICLETYPE>& contactDetection);
  ContactDetection<T, PARTICLETYPE>* getContactDetection();

  /// Particle-Fluid interaction for subgrid scale particles
  //  template<template<typename V> class DESCRIPTOR>
  //  void particleOnFluid(BlockLatticeStructure3D<T, DESCRIPTOR>& bLattice,
  //                       Cuboid3D<T>& cuboid, int overlap, T eps,
  //                       BlockGeometryStructure3D<T>& bGeometry);
  //  template<template<typename V> class DESCRIPTOR>
  //  void resetFluid(BlockLatticeStructure3D<T, DESCRIPTOR>& bLattice,
  //                  Cuboid3D<T>& cuboid, int overlap);

  friend class SuperParticleSystem3D<T, PARTICLETYPE>;
  friend class SuperParticleSysVtuWriter<T, PARTICLETYPE>;
  friend class SuperParticleSysVtuWriterMag<T>;
  friend class SimulateParticles<T, PARTICLETYPE>;

  //std::map<T, int> radiusDistribution();

  /// Integration method: explicit Euler
  void explicitEuler(T dT);

  ContactDetection<T, PARTICLETYPE>* getDetection()
  {
    return _contactDetection;
  }
  std::deque<PARTICLETYPE<T>> getParticles()
  {
    return _particles;
  }

protected:
  void integrateTorque(T dT);
  void integrateTorqueMag(T dT) {};
  void resetMag() {};

  void addShadowParticle(PARTICLETYPE<T>& p);

  mutable OstreamManager clout;
  SuperGeometry3D<T>& _superGeometry;
  LBconverter<T>& _converter;
  ContactDetection<T, PARTICLETYPE>* _contactDetection;
  SimulateParticles<T, PARTICLETYPE> _sim;

  std::deque<PARTICLETYPE<T> > _particles;
  std::deque<PARTICLETYPE<T> > _shadowParticles;
  std::list<std::shared_ptr<Force3D<T, PARTICLETYPE> > > _forces;
  std::list<std::shared_ptr<Boundary3D<T, PARTICLETYPE> > > _boundaries;

  std::vector<T> _physPos;
  std::vector<T> _physExtend;


  /// Integration methods, each need a special template particle
  void velocityVerlet1(T dT);
  void velocityVerlet2(T dT);
//  void implicitEuler(T dT, AnalyticalF3D<T,T>& getvel);
//  void adamBashforth4(T dT);
//  void predictorCorrector1(T dT);
//  void predictorCorrector2(T dT);
  void rungeKutta4_1(T dt);
  void rungeKutta4_2(T dt);
  void rungeKutta4_3(T dt);
  void rungeKutta4_4(T dt);
  void rungeKutta4(T dT);
  void updateParticleDistribution();
};

// Magnetic particle type
template<>
void ParticleSystem3D<double, MagneticParticle3D>::integrateTorqueMag(double dT);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::computeForce();
template<>
void ParticleSystem3D<double, MagneticParticle3D>::resetMag();

}  //namespace olb
#endif /* PARTICLE_SYSTEM_3D_H */
