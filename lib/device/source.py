import sys
import imp

particle = imp.load_source('particle', '/Users/aminmuha/Documents/myproject/cmi-injector/lib/particle/particle.py')

import numpy as np

class Source:
  def __init__(self, number_of_particles, muP, sigmaP, muV, sigmaV, distribution='gaussian'):
    self.particles = [] 
    self.number_of_particles = number_of_particles
    self.sigmaP = sigmaP
    self.sigmaV = sigmaV 
    self.muP = muP
    self.muV = muV
    self.distribution = distribution
    self.generate_particles()

  def generate_particles(self):
    position = np.random.normal(self.muP, self.sigmaP, self.number_of_particles)
    velocity = np.random.normal(self.muV, self.sigmaV, self.number_of_particles)
    for i in range(self.number_of_particles):
      self.particles.append(particle.SphericalParticle(0.00001, 0, 
           0, 0, 0, None, position=(position[i],0,0), velocity=(velocity[i],0,0)))



p = Source( 1000, 0, 3, 10, 5 )
for i in p.particles:
  print i.velocity, i. position
