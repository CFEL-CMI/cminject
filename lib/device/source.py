import sys
import imp

particle = imp.load_source('particle', '/Users/aminmuha/Documents/myproject/cmi-injector/lib/particle/particle.py')

import numpy as np

class Source:
  def __init__(self, number_of_particles, position, sigmaP, muV, sigmaV, distribution='gaussian'):
    self.particles = [] 
    self.number_of_particles = number_of_particles
    self.sigmaP = sigmaP
    self.sigmaV = sigmaV 
    self.muP = position
    self.muV = muV
    self.distribution = distribution
    self.generate_particles()

  def generate_particles(self):
    x = np.random.normal(self.muP[0], self.sigmaP[0], self.number_of_particles)
    y = np.random.normal(self.muP[1], self.sigmaP[1], self.number_of_particles) 
    z = np.random.normal(self.muP[2], self.sigmaP[2], self.number_of_particles)
    vx = np.random.normal(self.muV[0], self.sigmaV[0], self.number_of_particles)
    vy = np.random.normal(self.muV[1], self.sigmaV[1], self.number_of_particles)
    vz = np.random.normal(self.muV[2], self.sigmaV[2], self.number_of_particles)


    for i in range(self.number_of_particles):
      self.particles.append(particle.SphericalParticle(0.00001, 0, 
           0, 0, 0, None, position=(x[i], y[i]+0.01, z[i]), velocity=(vx[i],vy[i],vz[i])))



