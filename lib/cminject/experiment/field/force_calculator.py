from math import pi

class ForceCalculator:

  @staticmethod
  def DrageForceX(particle, fluid):
    force = 6 * pi * fluid.mu * particle.radius * (
		 fluid.fvx(particle.position) - particle.velocity[0] )
    return force

  @staticmethod
  def DrageForceY(particle, fluid):
    force = 6 * pi * fluid.mu * particle.radius * (
                 fluid.fvy(particle.position) - particle.velocity[1] )
    return force

  @staticmethod
  def DrageForceZ(particle, fluid):
    force = 6 * pi * fluid.mu * particle.radius * (
                 fluid.fvz(particle.position) - particle.velocity[2] )
    return force

  @staticmethod
  def PhotophoreticForce(particle, laser):
    return

  @staticmethod
  def RadiationPressure(particle, laser):
    return   
