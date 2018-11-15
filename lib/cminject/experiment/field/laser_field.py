from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad, dblquad
import numpy as np
from math import pi, sqrt, exp, cos, sin, tan, asin, acos, atan
import cmath

R = 8.314 # UNIVERSAL GAS CONSTANT

class Beam:
  def __init__(self, power, lamda, radius):
    self.power = power
    self.lamda = lamda
    self.w0 = radius

  def VortexIntensity( self, x, y, z):
    """Vortex intensity distribution propagating along the z axis"""
    z0 = pi * self.w0**2 / self.lamda    # Beyond the Rayleigh range of
    w = self.w0 * sqrt( 1 + (z/z0)**2 )
    r2 = x**2 + y**2
    I = ( 4 * self.power / pi ) * r2 / w**4 * exp( -2 * r2 / w**2 )
    return I


class CalculateAbsorption:
  def __init__(self, beam, particle):
    self.beam = beam
    self.nt = particle.index_of_ref
    self.position = particle.position
    self.r = particle.radius
    self.X=[]
    self.Y=[]
    self.Z=[]
    self.I=[]
    absorbed = dblquad(self.IntegratedAbsorbtion, 0, pi, lambda theta:   0, lambda theta:   2 * pi, args=('z'))
    self.absorbed = absorbed[0] * (self.r**2)

  def SplitSP(self, x, y, z, phi):
   """Split into s- and p-components. (x,y) position on the sphere
      phi angle between phase components of incoming beam"""

   r2 = x**2 + y**2 + z**2
   scomp = ((x * sin(phi) - y * cos(phi))**2) / r2
   pcomp = ((x * cos(phi) + y * sin(phi))**2) / r2
   return scomp, pcomp

  def AbsorptionCoeff(self, x, y):
   """calculates the absorption coefficients for S and P polarisation 
      components according to the Fresnel equations
      nt is the index of refraction of the sphere"""
   
   incident_angle = asin( sqrt(x**2 + y**2) )
   if incident_angle != 0:
     transmit_angle = cmath.asin( sin(incident_angle) / self.nt )
     delta_angle = incident_angle - transmit_angle
     total_angle = incident_angle + transmit_angle
     P = 1-abs( cmath.tan(delta_angle) / cmath.tan(total_angle))**2
     S = 1-abs( cmath.sin(delta_angle) / cmath.sin(total_angle))**2 
     return S, P, cos(incident_angle)
   if incident_angle == 0:
     P = 1 - abs((self.nt-1) / (self.nt+1) )**2
     return P, P, 1.0

  def PowerAbsorption(self, x, y, z):
    """ The absorption is calculated for each component s and p.
        the SplitSP function will return the percentage of S and P component in the Beam multiplied by beam intensity
        Then this should be multiplied by the absorption of S and P obtained by AbsiorptionCoeff function"""

    pabs, sabs, costheta = self.AbsorptionCoeff( x, y)
    phi = 0
    p, s = self.SplitSP(x, y, z, phi)
    pabs = pabs * p #* abs(z) 
    sabs = sabs * s #* abs(z)
    intensity = self.beam.VortexIntensity(self.position[0]+x, self.position[1]+y, self.position[2]+z)/self.r**2
    TotalAbsorpedPower = (intensity) * (pabs + sabs)
    return TotalAbsorpedPower, TotalAbsorpedPower*x, TotalAbsorpedPower*y, TotalAbsorpedPower*z

  def IntegratedAbsorbtion(self, theta, phi, comp): 
    """ integrate over a sphere"""
    x = sin(theta)*cos(phi) * self.r 
    y = sin(theta)*sin(phi) * self.r
    z = cos(theta)          * self.r
    self.X.append(x)
    self.Y.append(y)
    self.Z.append(z)
    
    power =  self.PowerAbsorption(x, y, z)
    if comp=='x':
       integ = abs(power[1]) * sin(phi)
    if comp=='z':
       integ = (power[3]) * sin(phi)
    self.I.append(integ/self.r**2)
    return integ


class CalculatePhotophoreticForces:
  def __init__(self, beam, particle, fluid):
    self.beam = beam
    self.particle = particle
    self.fluid = fluid
    D = (pi/2) * sqrt(pi/3) * fluid.k * fluid.eta * sqrt(fluid.T * 8 * R / (pi * fluid.mM)) / fluid.T
    CharP = (3/pi) * D * fluid.T / particle.radius
    first_half_sphere =  dblquad(self.IntegratedIntensity, 0, pi, lambda phi: 0, lambda phi:  pi)[0]
    second_half_sphere = dblquad(self.IntegratedIntensity, 0, pi, lambda phi: pi, lambda phi:  2*pi)[0]   
    first_half_sphere_ax =  dblquad(self.IntegratedIntensity, 0, pi/2, lambda phi: 0, lambda phi:  2*pi)[0]
    second_half_sphere_ax = dblquad(self.IntegratedIntensity, pi/2, pi, lambda phi: 0, lambda phi:  2*pi)[0]
    diff_trans = second_half_sphere - first_half_sphere
    diff_ax = second_half_sphere_ax - first_half_sphere_ax
    self.Fpp_tr = D * particle.radius**2 * fluid.P * diff_trans / (2 * particle.kp * CharP)
    self.Fpp_ax = D * particle.radius**2 * fluid.P * diff_ax / (2 * particle.kp * CharP)
#    print first_half_sphere, second_half_sphere, diff_trans, self.Fpp, self.Fpp_ax

 
  def IntegratedIntensity(self, theta, phi):
    x = sin(theta)*cos(phi) * self.particle.radius
    y = sin(theta)*sin(phi) * self.particle.radius
    z = cos(theta)          * self.particle.radius
    integ = sin(phi)*self.beam.VortexIntensity(self.particle.position[0]+x, self.particle.position[1]+y, self.particle.position[2]+z)
    return integ
