import sys
sys.path.insert(0, '../lib')
from scipy.integrate import ode
from field.interaction_field import *
from field.force_calculator import *
from als import AerodynamicsLensStack
from bgc import BufferGasCell

class Experiment:
  """This class carry out the info about the experiment you want to simulate. 
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment""" 

  def __init__(self, name, date, source, detector=None, devices=None, field=None, beam=None, traj=1000, end=0.8, dt=0.000001):
    self.name = name
    self.date = date
    self.field = field
    self.source = source
    self.devices = devices
    self.beam = beam
    self.traj = traj
    self.dt = dt
    self.NoTrajectories = (end/dt) / traj
    self.CalculateTrajectory(tStart=0, tEnd=end, dt=dt)
    self.SaveTrajectories()
 
  def CalculateTrajectory(self, tStart, tEnd, dt ):
    """ Calculate the particle trajectories by integrating the equation of motion"""
    count = 0
    for i in self.source.particles:
      traj = 0
      integral = ode( i.get_v_and_a )
      integral.set_integrator('lsoda') #, method='BDF',with_jacobian=False,atol=1e-6,rtol=1e-6,first_step=1e-5,nsteps=100000)
#      integral.set_integrator('vode',nsteps=500,method='bdf')
      integral.set_initial_value( (i.position + i.velocity), tStart ).set_f_params(self.field, self.beam)
      print "Calculate particle", count,  i.position
      while integral.successful() and self.ParticleInBoundary((integral.y[0], integral.y[1], integral.y[2]) ) and integral.t < tEnd:
        if traj%self.traj:
          i.trajectory.append( (integral.y[0], integral.y[1], integral.y[2]) )
          i.velocities.append( (integral.y[3], integral.y[4], integral.y[5]) )
        integral.integrate(integral.t + dt)
        i.position = (integral.y[0], integral.y[1], integral.y[2])
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
    if self.devices!=None:
      for i in self.devices:
        if type(i) is AerodynamicsLensStack:
          return i.ParticleInside(position)
        if type(i) is BufferGasCell:
          return i.BufferGasCellInside(position)
           
    else:
     #print "No Device Exists in The Experiment" 
     return True
  
  def SaveTrajectories(self):
    f = open(self.name+"_"+self.date,'w+')
    check = True
    for t in range(int(self.NoTrajectories)):
      if check:
        f.write(str(t*self.traj)+'\n')
        check = False # incase no particles propagating no need to write the time step
      c=0
      for i in self.source.particles:
        if len(i.trajectory)>t:  
          f.write(str(c) + " " + str(i.trajectory[t][0]) + " " + str(i.trajectory[t][1])
                      + " " + str(i.trajectory[t][2]) + " " + str(i.velocities[t][0]) + " "
                             + str(i.velocities[t][1]) + " " + str(i.velocities[t][2]) + '\n')
          c+=1
          check = True
