import cminject.experiment.particle as particle
import numpy as np

class Source:
  """This class generate particles with normal distribution of their initial position and momentum"""

  def __init__(self, number_of_particles, position, sigmaP, muV, sigmaV, radius=0.00001, rho=1, distribution='gaussian', seed=True):
    self.particles = [] #multiprocessing.Manager().list() 
    self.number_of_particles = number_of_particles
    self.sigmaP = sigmaP
    self.sigmaV = sigmaV 
    self.muP = position
    self.muV = muV
    self.distribution = distribution
    self.radius = radius
    self.rho = rho
    if seed:
      np.random.seed(1000)
    self.generate_particles()

  def generate_particles(self):
    """ Given the parameters of the normal distribution this function generates the intial positions
        and momentum of the particles"""

    if self.number_of_particles > 1:
      x = np.random.normal(self.muP[0], self.sigmaP[0], self.number_of_particles)
      y = np.random.normal(self.muP[1], self.sigmaP[1], self.number_of_particles) 
      z = np.random.normal(self.muP[2], self.sigmaP[2], self.number_of_particles)
      vx = np.random.normal(self.muV[0], self.sigmaV[0], self.number_of_particles)
      vy = np.random.normal(self.muV[1], self.sigmaV[1], self.number_of_particles)
      vz = np.random.normal(self.muV[2], self.sigmaV[2], self.number_of_particles)
    
    elif self.number_of_particles == 1:
      x = [self.muP[0]]
      y = [self.muP[1]]
      z = [self.muP[2]]
      vx = [self.muV[0]]
      vy = [self.muV[1]]
      vz = [self.muV[2]]

    elif self.number_of_particles == 0:
      raise ValueError('The number of particles cannot be zero')

    for i in range(self.number_of_particles): # create an object for each particle
      self.particles.append(particle.SphericalParticle(i, self.radius, self.rho, 
   			position=(x[i], y[i], z[i]), velocity=(vx[i],vy[i],vz[i])))



