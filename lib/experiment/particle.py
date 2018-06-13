import sys
sys.path.insert(0, '../lib')
from scipy.integrate import ode
from math import pi, atan, sin, cos
from field.laser_field import *
import numpy as np 

class Particle:
   def __init__(self, temp, mass, density, thermal_conductivity):
     self.temp = temp


class SphericalParticle(Particle):
  """This class define the characteristics of the particles that will be added to the simulation"""

  def __init__(self, radius, density, temp=298., thermal_conductivity=6.3, index_of_ref=0, 
                                             cp=1400., boundary=None, position=(0,0,0), velocity=(0,0,0)):
    self.radius = radius
    self.T = temp
    self.rho = density
    self.kp = thermal_conductivity     # 6.3 W/mK for carbon, 318 for Au. I made carbos as defult
    self.index_of_ref = index_of_ref
    self.position = position
    self.cp = cp # specific heat
    self.initialPosition = position
    self.velocity = velocity
    self.acceleration = (0, 0, 0)
    self.boundary = boundary
    self.M = self.mass()
    self.trajectory =[]
    self.velocities =[]
    self.insideFluid = True
    self.surfaceA = 4 * pi * radius**2
    self.collisions = 0
    self.called=0   # for debugging
    self.pressureAround = 0
    self.Tt = 298.
    self.time=0 

  def get_v_and_a(self, t, p_and_v, fluid, beam):
    """ This fuction returns the derivatives of the position and velocities for the integrator"""

    if fluid.FlowField:
      a = self.DragForceVector(fluid, p_and_v)/self.M
      self.CalculateTemp(fluid, t)
      self.CalculateCollisions(fluid, t)
      self.time+=t
#      print self.time, self.T, self.Tt
    if beam is not None:
      a[0] += self.FppTrans(fluid, beam)*cos(atan(self.position[1]/self.position[0])) / self.M
      a[1] += self.FppTrans(fluid, beam)*sin(atan(self.position[1]/self.position[0])) / self.M
      a[2] += self.FppAxial(fluid, beam) /  self.M
    return np.concatenate((p_and_v[3:], a))  # return the velocities and accelerations

  def DragForceVector(self,fluid, p_and_v):
     """This function calculates the drag force using Stokes' law for spherical particles in continuum"""
     try:
      token = fluid.fdrag((p_and_v[0],p_and_v[1],p_and_v[2]))
      fluid.vel = token[:3] - p_and_v[3:] # fuild velocity relative to the particle
      fluid.pressure = token[3]
      force_vector = 6 * pi * fluid.mu * self.radius * fluid.vel
      return force_vector/ self.SlipCorrection(fluid)
     except Exception as e:
       self.insideFluid = False  # The particle is outside the flow field
       return np.zeros(3)

  def SlipCorrection(self, fluid):
    """This function calculates the slip correction factor with temperature
       corrections. The Sutherland constant for helium is 79.4 at reference 
       Temperature of 273.K. I took a reference pressure of 1 Pascal, in 
       which the mean free path of helium is 0.01254. see J.Aersol Sci. 1976
       Vol. 7. pp 381-387 by Klaus Willeke"""

    Kn = 0.01254 * 1/(fluid.pressure*self.radius) * fluid.T/273. * (1 + 79.4/273.)/(1 + 79.4/fluid.T) 
    S = 1 + Kn * (1.246 + (0.42 * exp(-0.87/Kn)))
    return S


  def SlipCorrectionRoomTemp(self, fluid):
    """This function calculates the standard slip correction factor at room temperature"""
 
    lamda = fluid.MeanFreePath()
    Kn = lamda/self.radius
    S = 1 + Kn * (1.2310 + (0.4695 * exp(-1.1783/Kn)))
    return S

  def CalculateTemp(self, fluid, t):
    """This function calculates the particle temperature based on Newton Law of Cooling"""

    C = self.cp * self.M # The total heat capacity C of a system
    D = 2*self.radius
    muS = 1.96e-5        # dynamic viscosity at room temperature
    Re = fluid.Reynolds(D)
    Pr = fluid.Pr()
    M = fluid.Mach()
    Nu0 = fluid.Nusselt(D, muS)
    Nu = Nu0/ (1 + 3.42 * (M * Nu0/(Re * Pr)))
    h = fluid.h(Nu, D)
    r = h * self.surfaceA / C
    self.T = fluid.T + (self.T - fluid.T) * exp(-r*t)

  def CalculateCollisions(self, fluid, dt):
    """This function calculates the number of collisions based on Kinetic Gas Theory"""

    K = 1.38065e-23
    avG = 6.022e23
    v = sqrt(8*fluid.T*K/(fluid.mGas*pi))  # I should change this to velocity obtained from
    v += sqrt(fluid.vel[0]**2 + fluid.vel[1]**2 +fluid.vel[2]**2)
    n = fluid.pressure * avG/(8.314*fluid.T) # flow field
    C = 0.50 * n * v * self.surfaceA * dt
    k = (self.M + fluid.mGas)**2 / (2 * self.M * fluid.mGas)
    dT = (self.Tt - fluid.T)/k
    dT = C * dT
    self.Tt = fluid.T + (self.Tt - fluid.T) * exp(-C/k)
    self.collisions+=C
    
  def FppTrans(self, fluid, beam):
    """ This function calculates the trans component of the photophertic forces based on a semi imperical model"""

    if beam is not None and fluid is not None:
      fpp = CalculatePhotophoreticForces(beam, self, fluid)
      return fpp.Fpp_tr
    else:
      return 0

  def FppAxial(self, fluid, beam):
    """ This function calculates the axial component of the photophertic forces based on a semi imperical model"""
    if beam is not None and fluid is not None:
      fpp = CalculatePhotophoreticForces(beam, self, fluid)
      return fpp.Fpp_ax
    else:
      return 0

  def mass(self):
   """retuen the mass of a spherical particle"""

   return self.rho * 4/3 * pi * (self.radius)**3
