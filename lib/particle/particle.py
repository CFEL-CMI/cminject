class Particle:
   def __init__(self, temp, mass, density, thermal_conductivity):
     self.temp = temp


class SphericalParticle(Particle):
  """This class define the characteristics of the particles that will be added to the simulation"""
  def __init__(self, radius, temp, density, thermal_conductivity, index_of_ref ):
    self.radius = radius
    self.temp = temp
    self.density = density
    self.thermal_conductivity = thermal_conductivity 
    self.index_of_ref =index_of_ref
