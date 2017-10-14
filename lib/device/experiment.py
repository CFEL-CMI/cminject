import sys
sys.path.insert(0, '../lib')
from scipy.integrate import ode
from field.interaction_field import *
from field.force_calculator import *


class Experiment:
  """This class carry out the info about the experiment you want to simulate. 
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment""" 

  def __init__(self, name, date, source, detector=None, als=None, field=None):
    self.name = name
    self.date = date
    self.field = field
    self.device = []
    self.source = source
    self.als = als
    self.CalculateTrajectory()
    
  def CalculateTrajectory(self, tStart=0, tEnd=2 ):
    """ Calculate the particle trajectories by integrating the equation of motion"""
    count = 0
    for i in self.source.particles:
      traj = 0
      integral = ode( i.get_v_and_a )
      integral.set_integrator('lsoda')   #,method='BDF',with_jacobian=False,atol=1e-6,rtol=1e-6,first_step=1e-5,nsteps=100000)
      integral.set_initial_value( (i.position + i.velocity), tStart ).set_f_params(self.field)
      print "Calculate particle", count
      while integral.successful() and self.ParticleInBoundary((integral.y[0], integral.y[1], integral.y[2]) ) and integral.t < tEnd:
        if traj ==1:
          i.trajectory.append( (integral.y[0], integral.y[1], integral.y[2]) )
          traj = 0
        integral.integrate(integral.t + 0.00001)
        traj+=1
      count+=1
    return integral

  def LongestTrajectory(self):
    longest = 0
    for i in self.source.particles:
       if len(i.trajectory)>longest:
         longest = len(i.trajectory)
    return longest      

  def ParticleInBoundary(self, position):
    x, y, z = position
    inside = False
    offset = self.als.position[2]
    for i in self.als.segments:
       if (x**2 +y**2) < i[0]**2 and abs(z) >= abs(offset) and abs(z) < abs(offset+i[1]):
         return True
       offset = offset + i[1]
    return False
    
