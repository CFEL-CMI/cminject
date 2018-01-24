import sys
sys.path.insert(0, '../lib')
from scipy.integrate import ode
from field.interaction_field import *
from field.force_calculator import *
from als import AerodynamicsLensStack

class Experiment:
  """This class carry out the info about the experiment you want to simulate. 
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment""" 

  def __init__(self, name, date, source, detector=None, devices=None, field=None, beam=None):
    self.name = name
    self.date = date
    self.field = field
    self.source = source
    self.devices = devices
    self.beam = beam
    self.CalculateTrajectory()
    
  def CalculateTrajectory(self, tStart=0, tEnd=0.1 ):
    """ Calculate the particle trajectories by integrating the equation of motion"""
    count = 0
    for i in self.source.particles:
      traj = 0
      integral = ode( i.get_v_and_a )
      integral.set_integrator('lsoda')   #,method='BDF',with_jacobian=False,atol=1e-6,rtol=1e-6,first_step=1e-5,nsteps=100000)
      integral.set_initial_value( (i.position + i.velocity), tStart ).set_f_params(self.field, self.beam)
      print "Calculate particle", count
      while integral.successful() and self.ParticleInBoundary((integral.y[0], integral.y[1], integral.y[2]) ) and integral.t < tEnd:
#        if traj ==1:
        i.trajectory.append( (integral.y[0], integral.y[1], integral.y[2]) )
#          traj = 0
        integral.integrate(integral.t + 0.00001)
        i.position = (integral.y[0], integral.y[1], integral.y[2])
#        traj+=1
      count+=1
    return integral

  def LongestTrajectory(self):
    longest = 0
    for i in self.source.particles:
       if len(i.trajectory)>longest:
         longest = len(i.trajectory)
    return longest      

  def ParticleInBoundary(self, position):
    if self.devices!=None:
      for i in self.devices:
        if type(i) is AerodynamicsLensStack:
          x, y, z = position
          inside = False
          offset = i.position[2]
          for j in i.segments:
             if (x**2 +y**2) < j[0]**2 and abs(z) >= abs(offset) and abs(z) < abs(offset+j[1]):
               return True
             offset = offset + j[1]
          return False
      
    else:
     print "No Device Exists in The Experiment" 
     return True
    
