import time
from scipy.integrate import ode
from field.interaction_field import *
from field.force_calculator import *
from device.als import AerodynamicsLensStack
from device.bgc import BufferGasCell
from experiment.detector import Detector
from multiprocessing import Pool
from functools import partial
import numpy as np
import os


def SingeParticleTrajectory(i, count, field, beam, detector, tStart, tEnd, dt ):
  """ This function integrate the trajectory of a single particle. The propagation of particles in time stops 
      if the particle outside the boundary of the devices or the integration was not successful"""
  print(i)
  integral = ode( i.get_v_and_a )
  integral.set_integrator('lsoda', method='BDF',with_jacobian=False,atol=1e-8,rtol=1e-4,first_step=1e-5,nsteps=10000)
  integral.set_initial_value( (np.array(i.position + i.velocity)), tStart ).set_f_params(field, beam)
  print("Calculate particle", count,  i.position, i.velocity, i.ID)
  while integral.successful() and integral.t < tEnd and abs(integral.y[2]) < detector.end:
    integral.integrate(integral.t + dt)
    i.position = (integral.y[0], integral.y[1], integral.y[2])
    i.velocity = (integral.y[3], integral.y[4], integral.y[5])
  print("Final Position", i.position, i.velocity, '%2E' % i.collisions)



class Experiment:
  """This class carry out the info about the experiment you want to simulate.
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment"""

  def __init__(self, name, date, source, detector=None, devices=None, field=None, beam=None, traj=500, end=1.8, dt=1.e-5, directory='./', filename="v_and_p"):
    self.name = name
    self.date = date
    self.field = field
    self.source = source
    self.devices = devices
    self.beam = beam
    self.traj = traj
    self.dt = dt
    self.detector = detector
    self.NoTrajectories = (end/dt) / traj
    if not os.path.exists(directory):
      os.makedirs(directory)
    if not os.path.exists(directory+"/images/"):
      os.makedirs(directory+"/images/")

    self.directory = directory
    self.interp_time = 0
    self.called = 0
    
    self.CalculateTrajectory(tStart=0, tEnd=end, dt=dt)
#    self.SaveTrajectories()
#    self.detector.PlotDetector(directory, filename, source)

  def CalculateTrajectory(self, tStart, tEnd, dt ):
    """ Calculate the trajectories for all particles by integrating the equation of motion"""
    count = 1
#    for i in self.source.particles:
#      self.detector.ID = i.ID
#      SingeParticleTrajectory(i, count, self.field, self.beam, self.detector, tStart, tEnd, dt )
#      count+=1
    p = Pool()
    parallel_traj = partial(SingeParticleTrajectory, count=count, field=self.field, 
                         beam=self.beam, detector=self.detector, tStart=tStart, tEnd=tEnd, dt=dt)
    p.map(parallel_traj, self.source.particles)

    return True

  def IntegrateSingeParticle(self, i, count, tStart, tEnd, dt ):
    """ This function integrate the trajectory of a single particle. The propagation of particles in time stops 
        if the particle outside the boundary of the devices or the integration was not successful"""
    traj = 0
    integral = ode( i.get_v_and_a )
    integral.set_integrator('lsoda', method='BDF',with_jacobian=False,atol=1e-8,rtol=1e-4,first_step=1e-5,nsteps=10000)
    integral.set_initial_value( (np.array(i.position + i.velocity)), tStart ).set_f_params(self.field, self.beam)
    print("Calculate particle", count,  i.position, i.velocity, i.ID)
    while integral.successful() and self.ParticleInBoundary(i,integral.y[2]) and integral.t < tEnd and abs(integral.y[2]) < self.detector.end:
      if traj%self.traj:
        i.trajectory.append( (integral.y[0], integral.y[1], integral.y[2]) )
        i.velocities.append( (integral.y[3], integral.y[4], integral.y[5]) )
      self.CheckNearDetector(integral, i)
      integral.integrate(integral.t + dt)
      i.position = (integral.y[0], integral.y[1], integral.y[2])
      i.velocity = (integral.y[3], integral.y[4], integral.y[5])
      traj+=1
    try:
      self.detector.Projection()
    except:
      pass
    self.detector.x=[]; self.detector.y=[]; self.detector.z=[]
    self.detector.vx=[]; self.detector.vy=[]; self.detector.vz=[]
    self.detector.T=[]; self.detector.Tt=[]
    print("Final Position", i.position, i.velocity, '%2E' % i.collisions)


  def CheckNearDetector(self, integral, P):
    if self.detector.Check(abs(integral.y[2])):
        self.detector.T.append(P.T)
        self.detector.Tt.append(P.Tt)
        self.detector.x.append(integral.y[0])
        self.detector.y.append(integral.y[1])
        self.detector.z.append(integral.y[2])
        self.detector.vx.append(integral.y[3])
        self.detector.vy.append(integral.y[4])
        self.detector.vz.append(integral.y[5])
       

  """
  def FlyParticles(self, tStart, tEnd, dt):
    # This function integrates the trajecoties in parallel. It has faster interpolation but it is slower because the integrator
    #   uses a small time step for all particles. Thus, it is not used. I left it there in case I found a solution to this issue

    print("Calculate particle")
    integral = ode(self.get_v_a)
    integral.set_integrator('lsoda',  method='BDF',with_jacobian=False, atol=1e-6,rtol=1e-4,first_step=1e-5,nsteps=300)
    initial = np.concatenate((self.source.allParticlesPositions, self.source.allParticlesVelocities))
    integral.set_initial_value( initial.flatten(), tStart ) #.set_f_params()
    while integral.successful() and integral.t < tEnd:
        integral.integrate(integral.t + dt)
    print integral.y
    print 'interpolation: ',self.interp_time
    print "number of calls", self.called
    return integral

  def get_v_a(self,  t, p_and_v):
   self.called+=1
   p_and_v = p_and_v.reshape((2 * self.source.number_of_particles, 3))
   m = self.source.rho * 4/3 * pi * (self.source.radius)**3
   t1=time.time()
   a =  6 * pi * self.field.mu * self.source.radius * (
           self.field.fdrag(p_and_v[:self.source.number_of_particles])
                  - p_and_v[self.source.number_of_particles:] )/ (m * self.field.SlipCorrection)
   self.interp_time+=time.time()-t1
   v_a = np.concatenate((p_and_v[self.source.number_of_particles:], a))
   return v_a.flatten()
  """

  def LongestTrajectory(self):
    """ This function returns the longest trajectory of all particles and used only for trajectories visualizations"""

    longest = 0
    for i in self.source.particles:
       if len(i.trajectory)>longest:
         longest = len(i.trajectory)
    return longest

  def ParticleInBoundary(self, particle, Z):
#####################################
#    if self.devices==None:   # if there are no devices then use the boundary of the flow field
    if particle.insideFluid:
         return True
    else:
        if abs(Z)>=self.field.maxZ-self.field.maxZ/100.0:# or abs(Z)<=self.field.minZ: # if the particle is out but from the outlet it shouldn't stop
          return True               # till it reaches the detector
        else:
          return False
#####################################
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
    f = open(self.directory+self.name+"_"+self.date,'w+')
    check = True
    for t in range(int(self.NoTrajectories)):
      if check:
        f.write(str(t*self.traj)+'\n')
        check = False # incase no particles propagating no need to write the time step
      c=0
      for i in self.source.particles:
#        if len(i.trajectory)>t:
          f.write(str(c) + " " + str(i.trajectory[t][0]) + " " + str(i.trajectory[t][1])
                      + " " + str(i.trajectory[t][2]) + " " + str(i.velocities[t][0]) + " "
                             + str(i.velocities[t][1]) + " " + str(i.velocities[t][2]) + '\n')
          c+=1
          check = True
