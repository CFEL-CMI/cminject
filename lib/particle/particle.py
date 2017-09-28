
from scipy.integrate import ode

class Particle:
   def __init__(self, temp, mass, density, thermal_conductivity):
     self.temp = temp


class SphericalParticle(Particle):
  """This class define the characteristics of the particles that will be added to the simulation"""

  def __init__(self, radius, temp, density, thermal_conductivity, index_of_ref, boundary=None, position=(0,0,0), velocity=(0,0,0), acceleration=(0,0,0)):
    self.radius = radius
    self.temp = temp
    self.density = density
    self.thermal_conductivity = thermal_conductivity 
    self.index_of_ref = index_of_ref
    self.position = position
    self.velocity = velocity
    self.acceleration = acceleration
    self.boundary = boundary

  def CalculateTrajectory(self, tStart ):
    """ Calculate the particle trajectories by integrating the equation of motion"""

    integral = ode( velocity + acceleration )       
    integral.set_integrator('vode',method='BDF',with_jacobian=False,atol=1e-6,rtol=1e-6,first_step=1e-5,nsteps=100000)
    integral.set_initial_value( position + velocity, tStart )        
    while integral.successful() and self.InBoundary():
        integral.integrate(integral.t + 0.00001)    
    return integral

  def InBoundary(self):
    """ This function return true if the particle is still contained in the boundary of the problem"""
    if self.boundary == None:
       return True


  def CalculateDragForce(self, fluid):
    force=3*pi*mu*dParticle*velocityDifference/((rhoParticle*4/3*pi*(dParticle/2)**3)*correctionFactor)
    return force 
