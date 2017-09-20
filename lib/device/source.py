


class Source:
  def __init__(self, particle, number_of_particles, sigma, distribution='gaussian'):
    self.particle = particle 
    self.number_of_particles = number_of_particles
    self.sigma = sigma 
    self.distribution = distribution
