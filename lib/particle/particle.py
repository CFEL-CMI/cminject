
from scipy.integrate import ode

class Particle:
   def __init__(self, temp, mass, density, thermal_conductivity):
     self.temp = temp


class SphericalParticle(Particle):
  """This class define the characteristics of the particles that will be added to the simulation"""

  def __init__(self, radius, temp, density, thermal_conductivity, index_of_ref, boundary=None, position=(0,0,0), velocity=(0,0,0)):
    self.radius = radius
    self.temp = temp
    self.rho = density
    self.thermal_conductivity = thermal_conductivity 
    self.index_of_ref = index_of_ref
    self.position = position
    self.velocity = velocity
    self.acceleration = (0, 0, 0)
    self.boundary = boundary
    self.trajectory =[]

  def InBoundary(self):
    """ This function return true if the particle is still contained in the boundary of the problem"""
    if self.boundary == None:
       return True


  def get_v_and_a(self):
    return list(self.velocity + self.acceleration)


  def mass(self):
    return self.rho * 4/3 * pi * (self.radius)**3  
