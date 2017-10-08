from math import pi

class ForceCalculator:
  def __init__(self, particle):
    self.particle = particle

  @staticmethod
  def DrageForce(self, fluid):
    force = 6 * pi * fluid.mu * self.particle.radius * ( fluid.velocity - 
		self.particle.velocity ) / (( self.particle.rho * 4/3 * pi * 
				(self.particle.radius)**3) * correctionFactor)
    return force

  @staticmethod
  def PhotophoreticForce(self, laser):
    return

  @staticmethod
  def RadiationPressure(self, laser):
    return   
