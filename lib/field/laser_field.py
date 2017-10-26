from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad, dblquad
import numpy as np
from math import pi, sqrt, exp, cos, sin, tan, asin, acos, atan
import matplotlib.pyplot as plt

X=[]
Y=[]
I=[]


class LaserBeam:
  def __init__(self, power, lamda, radius, position):
    self.power = power
    self.lamda = lamda
    self.w0 = radius
    self.position = position
    self.absorbed = dblquad(self.IntegratedAbsorbtion, 0, pi, lambda theta:   0, lambda theta:   2*pi, args=((0,0,0), 2.6))
    print self.absorbed
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(-0.00003, 0.00003)
    ax.set_ylim(-0.00003, 0.00003)
    ax.scatter(X, Y, c=I, cmap='viridis')
    plt.show()

  def VortexIntensity( self, x, y, z):
    """Vortex intensity distribution propagating along the z axis"""
    z0 = pi * self.w0**2 / self.lamda    # Beyond the Rayleigh range of
    w = self.w0 * sqrt( 1 + (z/z0)**2 )
    r2 = x**2 + y**2
    I = ( 4 * self.power / pi ) * r2 / w**4 * exp( -2 * r2 / w**2 )
#    print I
    return I

  def SplitSP(self, x, y, phi):
   """Split into s- and p-components. (x,y) position on the sphere
      phi angle between phase components of incoming beam"""

   r2 = x**2 + y**2
   scomp = ((x * sin(phi) - y * cos(phi))**2) / r2
   pcomp = ((x * cos(phi) + y * sin(phi))**2) / r2
   return scomp, pcomp

  def AbsorptionCoeff(self, x, y, nt):
   """calculates the absorption coefficients for S and P polarisation 
      components according to the Fresnel equations
      nt is the index of refraction of the sphere"""
   
   incident_angle = asin( sqrt(x**2 + y**2) )
   if incident_angle != 0:
     transmit_angle = asin( sin(incident_angle) / nt )
     delta_angle = incident_angle - transmit_angle
     total_angle = incident_angle + transmit_angle
     P = 1-abs( tan(delta_angle) / tan(total_angle))**2
     S = 1-abs( sin(delta_angle) / sin(total_angle))**2 
     return S, P, cos(incident_angle)
   if incident_angle == 0:
     P = 1 - abs((nt-1) / (nt+1) )**2
     return p, p, 1.0

  def PowerAbsorption(self, x, y, z, particle, nt):
    """ The absorption is calculated for each component s and p.
        the SplitSP function will return the percentage of S and P component in the Beam multiplied by beam intensity
        Then this should be multiplied by the absorption of S and P obtained by AbsiorptionCoeff function"""

    pabs, sabs, costheta = self.AbsorptionCoeff( x, y, nt)
    phi = 0
    p, s = self.SplitSP(x, y, phi)
    pabs = pabs * p * costheta
    sabs = sabs * s * costheta
    intensity = self.VortexIntensity(particle[0]+x, particle[1]+y, particle[2]+z) / 10**4
    TotalAbsorpedPower = (intensity/(0.000035)**2) * (pabs + sabs)
    return TotalAbsorpedPower

  def IntegratedAbsorbtion(self, theta, phi, particle, nt): 
    """ integrate over a sphere"""
    x = 0.000035 * sin(theta)*cos(phi) 
    y = 0.000035 * sin(theta)*sin(phi) 
    z = 0.000035 * cos(theta)
    integ = (0.000035)**2 * self.PowerAbsorption(x, y, z, particle, nt) * sin(phi) 
    X.append(x)
    Y.append(y)
    I.append(integ)
    return integ



