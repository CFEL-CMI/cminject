/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn
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

#ifndef PARTICLESYSTEM_3D_HH
#define PARTICLESYSTEM_3D_HH

#include <list>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <memory>
#include <stdexcept>

#include "particle3D.h"
#include "particleSystem3D.h"
//#include "../functors/frameChangeF3D.h" //check
//#include "../utilities/vectorHelpers.h" //check

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
ParticleSystem3D<T, PARTICLETYPE>::ParticleSystem3D(
  SuperGeometry3D<T>& superGeometry, LBconverter<T>& conv)
  : clout(std::cout, "ParticleSystem3D"),
    _superGeometry(superGeometry),
    _converter(conv),
    _contactDetection(new ContactDetection<T, PARTICLETYPE>(*this)),
    _sim(this)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
ParticleSystem3D<T, PARTICLETYPE>::ParticleSystem3D(
  const ParticleSystem3D<T, PARTICLETYPE>& pS)
  : clout(std::cout, "ParticleSystem3D"),
    _superGeometry(pS._superGeometry),
    _converter(pS._converter),
    _contactDetection(new ContactDetection<T, PARTICLETYPE>(*this)),
    _sim(this),
    _physPos(pS._physPos),
    _physExtend(pS._physExtend)
{
  _particles = pS._particles;
  _forces = pS._forces;
}

template<typename T, template<typename U> class PARTICLETYPE>
ParticleSystem3D<T, PARTICLETYPE>::ParticleSystem3D(
  ParticleSystem3D<T, PARTICLETYPE> && pS)
  :     clout(std::cout, "ParticleSystem3D"),
        _superGeometry(pS._superGeometry),
        _converter(pS._converter),
        _contactDetection(new ContactDetection<T, PARTICLETYPE>(*this)),
        _sim(this),
        _physPos(pS._physPos),
        _physExtend(pS._physExtend)
{
  _particles = std::move(pS._particles);
  _forces = std::move(pS._forces);
}

template<typename T, template<typename U> class PARTICLETYPE>
ContactDetection<T, PARTICLETYPE>* ParticleSystem3D<T, PARTICLETYPE>::getContactDetection()
{
  return _contactDetection;
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::setContactDetection(ContactDetection<T, PARTICLETYPE>& contactDetection)
{
  delete _contactDetection;
  _contactDetection = contactDetection.generate(*this);
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::setPosExt(std::vector<T> physPos,
    std::vector<T> physExtend)
{
  _physPos = physPos;
  _physExtend = physExtend;
}

template<typename T, template<typename U> class PARTICLETYPE>
const std::vector<T>& ParticleSystem3D<T, PARTICLETYPE>::getPhysPos()
{
  return _physPos;
}

template<typename T, template<typename U> class PARTICLETYPE>
const std::vector<T>& ParticleSystem3D<T, PARTICLETYPE>::getPhysExtend()
{
  return _physExtend;
}

template<typename T, template<typename U> class PARTICLETYPE>
PARTICLETYPE<T>& ParticleSystem3D<T, PARTICLETYPE>::operator[](const int i)
{
  if (i < (int) _particles.size()) {
    return _particles[i];
  } else if (i < (int) (_particles.size() + _shadowParticles.size())) {
    return _shadowParticles[i - _particles.size()];
  } else {
    std::ostringstream os;
    os << i;
    std::string err = "Cannot access element: " + os.str()
                      + " to large";
    throw std::runtime_error(err);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
const PARTICLETYPE<T>& ParticleSystem3D<T, PARTICLETYPE>::operator[](
  const int i) const
{
  if ((unsigned)i < _particles.size()) {
    return _particles[i];
  } else if (i - _particles.size() < _shadowParticles.size()) {
    return _shadowParticles[i - _particles.size()];
  } else {
    std::ostringstream os;
    os << i;
    std::string err = "Cannot access element: " + os.str()
                      + " to large";
    throw std::runtime_error(err);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
int ParticleSystem3D<T, PARTICLETYPE>::size()
{
  return _particles.size();
}

template<typename T, template<typename U> class PARTICLETYPE>
int ParticleSystem3D<T, PARTICLETYPE>::sizeInclShadow() const
{
  return _particles.size() + _shadowParticles.size();
}

template<typename T, template<typename U> class PARTICLETYPE>
int ParticleSystem3D<T, PARTICLETYPE>::numOfActiveParticles()
{
  int activeParticles = 0;
  for (auto p : _particles) {
    if (p.getActive()) {
      activeParticles++;
    }
  }
  return activeParticles;
}
/*
template<typename T, template<typename U> class PARTICLETYPE>
std::map<T, int> ParticleSystem3D<T, PARTICLETYPE>::radiusDistribution()
{
  std::map<T, int> radDistr;
  typename std::map<T, int>::iterator it;
  T rad = 0;
  for (auto p : _particles) {
    if (!p.getAggl()) {
      rad = p.getRad();
      if (radDistr.count(rad) > 0) {
        it = radDistr.find(rad);
        it->second = it->second + 1;
      } else {
        radDistr.insert(std::pair<T, int>(rad, 1));
      }
    }
  }
  return radDistr;
}
*/
template<typename T, template<typename U> class PARTICLETYPE>
int ParticleSystem3D<T, PARTICLETYPE>::numOfForces()
{
  return _forces.size();
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::addParticle(PARTICLETYPE<T> &p)
{
  _particles.push_back(p);
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::clearParticles()
{
  _particles.clear();
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::addShadowParticle(PARTICLETYPE<T> &p)
{
  _shadowParticles.push_back(p);
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::addForce(
  std::shared_ptr<Force3D<T, PARTICLETYPE> > pF)
{
  _forces.push_back(pF);
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::addBoundary(
  std::shared_ptr<Boundary3D<T, PARTICLETYPE> > b)
{
  _boundaries.push_back(b);
}

template<typename T, template<typename U> class PARTICLETYPE>
template<template<typename V> class DESCRIPTOR>
void ParticleSystem3D<T, PARTICLETYPE>::setVelToFluidVel(
  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> & fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::setVelToAnalyticalVel(
  AnalyticalConst3D<T,T>& aVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      std::vector<T> vel(3, T());
      aVel(&p.getVel()[0], &p.getPos()[0]);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::computeForce()
{
  typename std::deque<PARTICLETYPE<T> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
    if (p->getActive()) {
      p->resetForce();
      for (auto f : _forces) {
        f->applyForce(p, pInt, *this);
      }
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::computeBoundary()
{
  typename std::deque<PARTICLETYPE<T> >::iterator p;
  for (auto f : _boundaries) {
    for (p = _particles.begin(); p != _particles.end(); ++p) {
      if (p->getActive()) {
        f->applyBoundary(p, *this);
      }
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::simulate(T dT)
{
  _sim.simulate(dT);
}

template<typename T, template<typename U> class PARTICLETYPE>
int ParticleSystem3D<T, PARTICLETYPE>::countMaterial(int mat)
{
  int num = 0;
  int locLatCoords[3] = {0};
  for (auto& p : _particles) {
    _superGeometry.getCuboidGeometry().get(p.getCuboid()).getFloorLatticeR(
      locLatCoords, &p.getPos()[0]);
    const BlockGeometryStructure3D<T>& bg = _superGeometry.getExtendedBlockGeometry(_superGeometry.getLoadBalancer().loc(p.getCuboid()));
    int iX = locLatCoords[0]+_superGeometry.getOverlap();
    int iY = locLatCoords[1]+_superGeometry.getOverlap();
    int iZ = locLatCoords[2]+_superGeometry.getOverlap();
    if    (bg.get(iX,     iY,     iZ) == mat
           || bg.get(iX,     iY + 1, iZ) == mat
           || bg.get(iX,     iY,     iZ + 1) == mat
           || bg.get(iX,     iY + 1, iZ + 1) == mat
           || bg.get(iX + 1, iY,     iZ) == mat
           || bg.get(iX + 1, iY + 1, iZ) == mat
           || bg.get(iX + 1, iY,     iZ + 1) == mat
           || bg.get(iX + 1, iY + 1, locLatCoords[2] + 1) == mat) {
      num++;
    }
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::explicitEuler(T dT)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      for (int i = 0; i < 3; i++) {
        p.getVel()[i] += p.getForce()[i] * p.getInvMass() * dT;
        p.getPos()[i] += p.getVel()[i] * dT;

        // if particles are too fast, e.g. material boundary can not work anymore
#ifdef OLB_DEBUG
        if (p.getVel()[i] * dT > _superGeometry.getCuboidGeometry().getMaxDeltaR()) {
          std::cout << " PROBLEM: particle speed too high rel. to delta of "
                    "lattice: "<< std::endl;
          std::cout << "p.getVel()[i]*dT: " << i <<" "<< p.getVel()[i] * dT;
          std::cout << "MaxDeltaR(): " <<
                    _superGeometry.getCuboidGeometry().getMaxDeltaR() << std::endl;
          exit(-1);
        }
#endif
      }
    }
  }
}


//template<typename T, template<typename U> class PARTICLETYPE>
//void ParticleSystem3D<T, PARTICLETYPE>::integrateTorque(T dT)
//{
//  for (auto& p : _particles) {
//    if (p.getActive()) {
//      for (int i = 0; i < 3; i++) {
//        p.getAVel()[i] += p.getTorque()[i] * 1.
//                          / (2. / 5. * p.getMass() * std::pow(p.getRad(), 2)) * dT;
//      }
//    }
//  }
//}

//template<typename T, template<typename U> class PARTICLETYPE>
//void ParticleSystem3D<T, PARTICLETYPE>::integrateTorqueMag(T dT) {
////template<typename T>
////void ParticleSystem3D<T, MagneticParticle3D>::integrateTorqueMag(T dT) {
//  for (auto& p : _particles) {
//    if (p.getActive()) {
//      Vector<T, 3> deltaAngle;
//      T angle;
//      T epsilon = std::numeric_limits<T>::epsilon();
//      for (int i = 0; i < 3; i++) {
//        // change orientation due to torque moments
//        deltaAngle[i] = (5. * p.getTorque()[i] * dT * dT) / (2. * p.getMass() * std::pow(p.getRad(), 2));
//        // apply change in angle to dMoment vector
//      }
//      angle = norm(deltaAngle);
//      if (angle > epsilon) {
//        //Vector<T, 3> axis(deltaAngle);
//        //  Vector<T, 3> axis(T(), T(), T(1));
//        //axis.normalize();
//        std::vector<T> null(3, T());
//
//        //RotationRoundAxis3D<T, S> rotRAxis(p.getPos(), fromVector3(axis), angle);
//        RotationRoundAxis3D<T, S> rotRAxis(null, fromVector3(deltaAngle), angle);
//        T input[3] = {p.getMoment()[0], p.getMoment()[1], p.getMoment()[2]};
//        Vector<T, 3> in(input);
//        T output[3] = {T(), T(), T()};
//        rotRAxis(output, input);
//        Vector<T, 3> out(output);
////        std::vector<T> mainRefVec(3, T());
////        mainRefVec[0] = 1.;
////        std::vector<T> secondaryRefVec(3, T());
////        secondaryRefVec[2] = 1.;
////        AngleBetweenVectors3D<T, T> checkAngle(mainRefVec, secondaryRefVec);
////        T angle[1];
////        checkAngle(angle, input);
//        std::cout<< "|moment|_in: " << in.norm() << ", |moment|_out: " << out.norm()
//                 << ", | |in| - |out| |: " << fabs(in.norm() - out.norm())
//                 /*<< " Angle: " << angle[0] */<< std::endl;
//        p.getMoment()[0] = output[0];
//        p.getMoment()[1] = output[1];
//        p.getMoment()[2] = output[2];
//      }
//    }
//  }
//}

//template<typename T, template<typename U> class PARTICLETYPE>
//void ParticleSystem3D<T, PARTICLETYPE>::implicitEuler(T dT, AnalyticalF3D<T,T>& getvel) {
//  _activeParticles = 0;
//  for (auto& p : _particles) {
//    if(p.getActive()) {
//      std::vector<T> fVel = getvel(p._pos);
////      std::vector<T> vel = p.getVel();
//      std::vector<T> pos = p.getPos();
//      T C = 6.* M_PI * p._rad * this->_converter.getCharNu() * this->_converter.getCharRho()* dT/p.getMass();
//      for (int i = 0; i<3; i++) {
//        p._vel[i] = (p._vel[i]+C*fVel[i]) / (1+C);
//        p._pos[i] += p._vel[i] * dT;
//      }
////      p.setVel(vel);
////      p.setPos(pos);
//      checkActive(p);
////      cout << "C: " << C << std::endl;
//    }
////    cout << "pos: " << p._pos[0] << " " << p._pos[1] << " " << p._pos[2] << " " << p._vel[0] << " " << p._vel[1] << " " << p._vel[2]<< std::endl; //"\n force: " << p._force[0] << " " << p._force[1] << " " << p._force[2]<< "\n " << std::endl;
//  }
//}

/*
template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::predictorCorrector1(T dT)
{
  std::vector<T> vel;
  std::vector<T> pos;
  std::vector<T> frc;
  for (auto& p : _particles) {
    if (p.getActive()) {
      vel = p.getVel();
      p.setVel(vel, 1);
      pos = p.getPos();
      p.setVel(pos, 2);
      frc = p.getForce();
      p.setForce(frc, 1);
      for (int i = 0; i < 3; i++) {
        vel[i] += p._force[i] / p.getMass() * dT;
        pos[i] += vel[i] * dT;
      }
      p.setVel(vel);
      p.setPos(pos);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::predictorCorrector2(T dT)
{
  std::vector<T> vel;
  std::vector<T> pos;
  std::vector<T> frc;
  for (auto& p : _particles) {
    if (p.getActive()) {
      vel = p.getVel(1);
      pos = p.getVel(2);
      for (int i = 0; i < 3; i++) {
        vel[i] += dT * .5 * (p.getForce()[i] + p.getForce(1)[i]) / p.getMass();
        pos[i] += vel[i] * dT;
      }
      p.setVel(vel);
      p.setPos(pos);
    }
  }
}
*/

/*
template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::adamBashforth4(T dT)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      std::vector<T> vel = p.getVel();
      std::vector<T> pos = p.getPos();
      for (int i = 0; i < 3; i++) {
        vel[i] += dT / p.getMas()
                  * (55. / 24. * p.getForce()[i] - 59. / 24. * p.getForce(1)[i]
                     + 37. / 24. * p.getForce(2)[i] - 9. / 24. * p.getForce(3)[i]);
      }
      p.rotAndSetVel(vel);
      for (int i = 0; i < 3; i++) {
        pos[i] += dT
                  * (55. / 24. * p.getVel()[i] - 59. / 24. * p.getVel(1)[i]
                     + 37. / 24. * p.getVel(2)[i] - 3. / 8. * p.getVel(3)[i]);
      }
      p.setPos(pos);
    }
  }
}
*/

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::velocityVerlet1(T dT)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      std::vector<T> pos = p.getPos();
      for (int i = 0; i < 3; i++) {
        pos[i] += p.getVel()[i] * dT
                  + .5 * p.getForce()[i] / p.getMass() * std::pow(dT, 2);
      }
      p.setPos(pos);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::velocityVerlet2(T dT)
{
  std::vector<T> frc;
  std::vector<T> vel;
  typename std::deque<PARTICLETYPE<T> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
//  for (auto& p : _particles) {
    if (p->getActive()) {
      frc = p->getForce();
      vel = p->getVel();
      p->resetForce();
      for (auto f : _forces) {
        f->applyForce(p, pInt, *this);
      }
      for (int i = 0; i < 3; i++) {
        vel[i] += .5 * (p->getForce()[i] + frc[i]) / p->getMass() * dT;
      }
      p->setVel(vel);
    }
  }
}

/*
template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::rungeKutta4(T dT)
{
  rungeKutta4_1(dT);
  rungeKutta4_2(dT);
  rungeKutta4_3(dT);
  rungeKutta4_4(dT);
}


template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::rungeKutta4_1(T dT)
{
  typename std::deque<PARTICLETYPE<T> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
    if (p->getActive()) {
      p->resetForce();
      for (auto f : _forces) {
        f->applyForce(p, pInt, *this);
      }

      std::vector<T> k1 = p->getForce();
      p->setForce(k1, 3);
      std::vector<T> storeVel = p->getVel();
      p->setVel(storeVel, 3);
      std::vector<T> storePos = p->getPos();
      p->setVel(storePos, 2);

      std::vector<T> vel(3, T()), pos(3, T());
      for (int i = 0; i < 3; i++) {
        vel[i] = p->getVel(3)[i] + dT / 2. / p->getMass() * k1[i];
        pos[i] = p->getVel(2)[i] + dT / 2. * vel[i];
      }
      p->setVel(vel);
      p->setPos(pos);
    }
  }
}


template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::rungeKutta4_2(T dT)
{
  typename std::deque<PARTICLETYPE<T> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
    if (p->getActive()) {
      p->resetForce();
      for (auto f : _forces) {
        f->applyForce(p, pInt, *this);
      }
      std::vector<T> k2 = p->getForce();
      p->setForce(k2, 2);

      std::vector<T> vel(3, T()), pos(3, T());
      for (int i = 0; i < 3; i++) {
        vel[i] = p->getVel(3)[i] + dT / 2. / p->getMass() * k2[i];
        pos[i] = p->getVel(2)[i] + dT / 2. * vel[i];
      }
      p->setVel(vel);
      p->setPos(pos);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::rungeKutta4_3(T dT)
{
  typename std::deque<PARTICLETYPE<T> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
    if (p->getActive()) {
      p->resetForce();
      for (auto f : _forces) {
        f->applyForce(p, pInt, *this);
      }
      std::vector<T> k3 = p->getForce();
      p->setForce(k3, 1);
      std::vector<T> vel(3, T()), pos(3, T());
      for (int i = 0; i < 3; i++) {
        vel[i] = p->getVel(3)[i] + dT / p->getMass() * k3[i];
        pos[i] = p->getVel(2)[i] + dT * vel[i];
      }
      p->setVel(vel);
      p->setPos(pos);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::rungeKutta4_4(T dT)
{
  typename std::deque<PARTICLETYPE<T> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
    if (p->getActive()) {
      p->resetForce();
      for (auto f : _forces) {
        f->applyForce(p, pInt, *this);
      }
      std::vector<T> k4 = p->getForce();

      std::vector<T> vel(3, T()), pos(3, T());
      for (int i = 0; i < 3; i++) {
        vel[i] = p->getVel(3)[i]
                 + dT / 6. / p->getMass()
                 * (p->getForce(3)[i] + 2 * p->getForce(2)[i]
                    + 2 * p->getForce(1)[i] + p->getForce()[i]);
        pos[i] = p->getVel(2)[i] + dT * vel[i];
      }
      p->setVel(vel);
      p->setPos(pos);
    }
  }
}
*/


template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::saveToFile(std::string fullName)
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  std::vector<T> serPar(PARTICLETYPE<T>::serialPartSize, 0);
  for (auto& p : _particles) {
    p.serialize(&serPar[0]);
    typename std::vector<T>::iterator it = serPar.begin();
    for (; it != serPar.end(); ++it) {
      fout << *it << " ";
    }
    fout << std::endl;
    //fout << p.getPos()[0] << " " << p.getPos()[1] << " " << p.getPos()[2] << " " << p.getVel()[0] << " " << p.getVel()[1] << " "<< p.getVel()[2] << std::endl;
  }
  fout.close();
}

//template<typename T, template<typename U> class PARTICLETYPE>
//template<template<typename V> class DESCRIPTOR>
//void ParticleSystem3D<T, PARTICLETYPE>::particleOnFluid(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& bLattice, Cuboid3D<T>& cuboid,
//  int overlap, T eps, BlockGeometryStructure3D<T>& bGeometry)
//{
//  T rad = 0;
//  std::vector<T> vel(3, T());
//  T minT[3] = {0}, maxT[3] = {0}, physR[3] = {0};
//  int min[3] = {0}, max[3] = {0};
//  //cout << "OVERLAP: " << overlap << std::endl;
//  for (int i = 0; i < sizeInclShadow(); ++i) {
//    rad = operator[](i).getRad();
//    minT[0] = operator[](i).getPos()[0] - rad * eps;
//    minT[1] = operator[](i).getPos()[1] - rad * eps;
//    minT[2] = operator[](i).getPos()[2] - rad * eps;
//    maxT[0] = operator[](i).getPos()[0] + rad * eps;
//    maxT[1] = operator[](i).getPos()[1] + rad * eps;
//    maxT[2] = operator[](i).getPos()[2] + rad * eps;
//    cuboid.getLatticeR(min, minT);
//    cuboid.getLatticeR(max, maxT);
//    max[0]++;
//    max[1]++;
//    max[2]++;
//    T porosity = 0;
//    T dist = 0;
//    for (int iX = min[0]; iX < max[0]; ++iX) {
//      for (int iY = min[1]; iY < max[1]; ++iY) {
//        for (int iZ = min[2]; iZ < max[2]; ++iZ) {
//          cuboid.getPhysR(physR, iX, iY, iZ);
//          dist = std::pow(physR[0] - operator[](i).getPos()[0], 2)
//                 + std::pow(physR[1] - operator[](i).getPos()[1], 2)
//                 + std::pow(physR[2] - operator[](i).getPos()[2], 2);
//          if (dist < rad * rad * eps * eps) {
//            vel = operator[](i).getVel();
//            vel[0] = _converter.latticeVelocity(vel[0]);
//            vel[1] = _converter.latticeVelocity(vel[1]);
//            vel[2] = _converter.latticeVelocity(vel[2]);
//            if (dist < rad * rad) {
//              porosity = 0;
//              bLattice.get(iX, iY, iZ).defineExternalField(
//                DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, &porosity);
//              bLattice.get(iX, iY, iZ).defineExternalField(
//                DESCRIPTOR<T>::ExternalField::localDragBeginsAt,
//                DESCRIPTOR<T>::d, &vel[0]);
//            } else {
//              T d = std::sqrt(dist) - rad;
//              porosity = 1.
//                         - std::pow(std::cos(M_PI * d / (2. * (eps - 1.) * rad)), 2);
//              bLattice.get(iX, iY, iZ).defineExternalField(
//                DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, &porosity);
//              bLattice.get(iX, iY, iZ).defineExternalField(
//                DESCRIPTOR<T>::ExternalField::localDragBeginsAt,
//                DESCRIPTOR<T>::d, &vel[0]);
//            }
//          }
//        }
//      }
//    }
//  }
//}
//
//template<typename T, template<typename U> class PARTICLETYPE>
//template<template<typename V> class DESCRIPTOR>
//void ParticleSystem3D<T, PARTICLETYPE>::resetFluid(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& bLattice, Cuboid3D<T>& cuboid,
//  int overlap)
//{
//  T rad = 0, vel[3] = {0};
//  T minT[3] = {0}, maxT[3] = {0}, physR[3] = {0};
//  int min[3] = {0}, max[3] = {0};
//  for (int i = 0; i < sizeInclShadow(); ++i) {
//    rad = operator[](i).getRad();
//    minT[0] = operator[](i).getPos()[0] - rad * 1.2;
//    minT[1] = operator[](i).getPos()[1] - rad * 1.2;
//    minT[2] = operator[](i).getPos()[2] - rad * 1.2;
//    maxT[0] = operator[](i).getPos()[0] + rad * 1.2;
//    maxT[1] = operator[](i).getPos()[1] + rad * 1.2;
//    maxT[2] = operator[](i).getPos()[2] + rad * 1.2;
//    cuboid.getLatticeR(min, minT);
//    cuboid.getLatticeR(max, maxT);
//    max[0]++;
//    max[1]++;
//    max[2]++;
//    T porosity = 1;
//    T dist = 0;
//    for (int iX = min[0]; iX < max[0]; ++iX) {
//      for (int iY = min[1]; iY < max[1]; ++iY) {
//        for (int iZ = min[2]; iZ < max[2]; ++iZ) {
//          bLattice.get(iX, iY, iZ).defineExternalField(
//            DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, &porosity);
//          bLattice.get(iX, iY, iZ).defineExternalField(
//            DESCRIPTOR<T>::ExternalField::localDragBeginsAt, DESCRIPTOR<T>::d,
//            vel);
//          //          }
//        }
//      }
//    }
//  }
//}


template<>
void ParticleSystem3D<double, MagneticParticle3D>::resetMag()
{
  typename std::deque<MagneticParticle3D<double> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
    if (p->getActive()) {
      p->resetForce();
      p->resetTorque();
    }
  }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::computeForce()
{
  typename std::deque<MagneticParticle3D<double> >::iterator p;
  int pInt = 0;
  for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
    if (p->getActive()) {
      for (auto f : _forces) {
        f->applyForce(p, pInt, *this);
      }
    }
  }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::integrateTorqueMag(double dT)
{
  for (auto& p : _particles) {
    Vector<double, 3> deltaAngle;
    double angle;
    double epsilon = std::numeric_limits<double>::epsilon();
    double damping = std::pow((1. - p.getADamping()), dT);
    for (int i = 0; i < 3; i++) {
      p.getAVel()[i] += (5. * (p.getTorque()[i]) * dT) / (2.  * p.getMass() * std::pow(p.getRad(), 2));
      p.getAVel()[i] *= damping;
      deltaAngle[i] = p.getAVel()[i] * dT;
    }
    angle = norm(deltaAngle);
    if (angle > epsilon) {
      std::vector<double> null(3, double());

      RotationRoundAxis3D<double, double> rotRAxis(null, util::fromVector3(deltaAngle), angle);
      double input[3] = {p.getMoment()[0], p.getMoment()[1], p.getMoment()[2]};
      Vector<double, 3> in(input);
      double output[3] = {double(), double(), double()};
      rotRAxis(output, input);
      Vector<double, 3> out(output);
      // renormalize output
      if (out.norm() > epsilon) {
        out = (1. / out.norm()) * out;
      }

      p.getMoment()[0] = out[0];
      p.getMoment()[1] = out[1];
      p.getMoment()[2] = out[2];
    }
  }
}


}  //namespace olb
#endif /* PARTICLESYSTEM_3D_HH */
