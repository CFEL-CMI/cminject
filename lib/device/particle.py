import sys
sys.path.insert(0, '../lib')
from scipy.integrate import ode
from math import pi, atan, sin, cos
from field.laser_field import *

class Particle:
   def __init__(self, temp, mass, density, thermal_conductivity):
     self.temp = temp


class SphericalParticle(Particle):
  """This class define the characteristics of the particles that will be added to the simulation"""

  def __init__(self, radius, density, temp=0, thermal_conductivity=6.3, index_of_ref=0, boundary=None, position=(0,0,0), velocity=(0,0,0)):
    self.radius = radius
    self.temp = temp
    self.rho = density
    self.kp = thermal_conductivity     # 6.3 W/mK for carbon, 318 for Au. I made carbos as defult
    self.index_of_ref = index_of_ref
    self.position = position
    self.velocity = velocity
    self.acceleration = (0, 0, 0)
    self.boundary = boundary
    self.trajectory =[]
    self.velocities =[]

  def InBoundary(self):
    """ This function return true if the particle is still contained in the boundary of the problem"""
    if self.boundary == None:
       return True


  def get_v_and_a(self, t, p_and_v, fluid, beam):
    a=[]
    Fx = Fy = Fz = 0
    if fluid.FlowField:
      Fx = self.DrageForceX(fluid, p_and_v)
      Fy = self.DrageForceY(fluid, p_and_v)
      Fz = self.DrageForceZ(fluid, p_and_v)
    if beam is not None:
      Fx += self.FppTrans(fluid, beam)*cos(atan(self.position[1]/self.position[0]))
      Fy += self.FppTrans(fluid, beam)*sin(atan(self.position[1]/self.position[0]))
      Fz += self.FppAxial(fluid, beam)
#    print Fx, Fy, Fz
    a.append( Fx/self.mass())
    a.append( Fy/self.mass())
    a.append( Fz/self.mass())

    return list(list(p_and_v[3:]) + a)

  def DrageForceX(self, fluid, p_and_v):
    try:
      force = 6 * pi * fluid.mu * self.radius * (
                 fluid.fvx(p_and_v[:3]) - p_and_v[3] )
      return force
    except:
      return 0

  def DrageForceY(self, fluid, p_and_v):
    try:
      force = 6 * pi * fluid.mu * self.radius * (
                 fluid.fvy(p_and_v[:3]) - p_and_v[4] )
      return force
    except:
      return 0

  def DrageForceZ(self, fluid, p_and_v):
    try:
      force = 6 * pi * fluid.mu * self.radius * (
                 fluid.fvz(p_and_v[:3]) - p_and_v[5] )
      return force
    except:
      return 0

  def FppTrans(self, fluid, beam):
    if beam is not None and fluid is not None:
      fpp = CalculatePhotophoreticForces(beam, self, fluid)
      return fpp.Fpp_tr
    else:
      return 0
  def FppAxial(self, fluid, beam):
    if beam is not None and fluid is not None:
      fpp = CalculatePhotophoreticForces(beam, self, fluid) 
      return fpp.Fpp_ax
    else:
      return 0


  def mass(self):
    return self.rho * 4/3 * pi * (self.radius)**3  
