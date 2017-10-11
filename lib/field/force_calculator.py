from math import pi

class ForceCalculator:
  def __init__(self, particle):
    self.particle = particle

  @staticmethod
  def DrageForceX(self, fluid):
    force = 6 * pi * fluid.mu * self.particle.radius * (
		 fluid.fvx(self.particle.position) - self.particle.velocity )
    return force

  @staticmethod
  def DrageForceY(self, fluid):
    force = 6 * pi * fluid.mu * self.particle.radius * (
                 fluid.fvy(self.particle.position) - self.particle.velocity )
    return force

  @staticmethod
  def DrageForceZ(self, fluid):
    force = 6 * pi * fluid.mu * self.particle.radius * (
                 fluid.fvz(self.particle.position) - self.particle.velocity )
    return force

  @staticmethod
  def PhotophoreticForce(self, laser):
    return

  @staticmethod
  def RadiationPressure(self, laser):
    return   
