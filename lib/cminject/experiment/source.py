import cminject.experiment.particle as particle
import numpy as np

class Source:
  """This class generate particles with normal distribution of their initial position and momentum"""

  def __init__(self, number_of_particles, position, sigmaP, muV, sigmaV, 
               radius=1.e-6, sigma_radius=1.e-7, rho=1, distribution='gaussian', seed=None):
    self.particles = [] #multiprocessing.Manager().list() 
    self.number_of_particles = number_of_particles
    self.sigmaP = sigmaP
    self.sigmaV = sigmaV 
    self.muP = position
    self.muV = muV
    self.distribution = distribution
    self.radius = radius
    self.sigma_radius = sigma_radius
    self.rho = rho
    if seed is not None:
      np.random.seed(seed)
    self.generate_particles()

  def generate_particles(self):
    """ Given the parameters of the normal distribution this function generates the intial positions
        and momentum of the particles"""

    if self.number_of_particles >= 1:
      position = np.random.normal(self.muP, self.sigmaP,
                                  (self.number_of_particles, 3))
      velocity = np.random.normal(self.muV, self.sigmaV,
                                  (self.number_of_particles, 3))
      r = np.random.normal(self.radius, self.sigma_radius, self.number_of_particles)
    else:
      raise ValueError('The number of particles must be >= 1')

    for i in range(self.number_of_particles): # create an object for each particle
      self.particles.append(
        particle.SphericalParticle(
          i, r[i], self.rho,
          position=tuple(position[i]),
          velocity=tuple(velocity[i])
        )
      )



