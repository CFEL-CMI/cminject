import traceback
from scipy.integrate import ode
from multiprocessing import Pool, Manager, cpu_count
from functools import partial
import numpy as np
import h5py
import os
import sys

def SingeParticleTrajectory(particle, devices, field, beam, detector, tStart, tEnd, dt, zlim, l, traj=False):
  """ This function integrate the trajectory of a single particle. The propagation of particles in time stops 
      if the particle outside the boundary of the devices or the integration was not successful"""
  integral = ode( particle.get_v_and_a )
  integral.set_integrator('lsoda') #, method='BDF',with_jacobian=False,atol=1e-8,rtol=1e-4,first_step=1e-5,nsteps=10000)
  integral.set_initial_value( np.concatenate((particle.position, 
                              particle.velocity)), tStart ).set_f_params(field, beam)

  if traj:    # If trajectory is true then store the initial position
    particle.trajectory.append([0., particle.position[0], particle.position[1], particle.position[2], 
                                particle.velocity[0], particle.velocity[1], particle.velocity[2], 
                                particle.T, particle.Tt])  

  while integral.successful() and integral.t < tEnd:
    if particle.CheckParticleIn(integral.y, zlim, detector.end, devices)==0:
      if traj:
        particle.trajectory.append([integral.t, integral.y[0], integral.y[1], integral.y[2],
                             integral.y[3], integral.y[4], integral.y[5], particle.T, particle.Tt])
      if  field is not None and field.pressure>0:
        particle.CalculateTemp(field, dt)
        particle.CalculateCollisions(field, dt)
        if particle.T > 77.:
            particle.TimeLN[0] = integral.t
        if particle.Tt > 77.:
            particle.TimeLN[1] = integral.t
      particle.position = integral.y[:3]
      particle.velocity = integral.y[3:]
      particle.TOF = integral.t     
      integral.integrate(integral.t + dt)
    else:
      break
  l.append(particle)

class Experiment:
  """This class carry the info about the experiment you want to simulate.
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment"""

  def __init__(self, name, date, source, detector=None, devices=None, field=None, 
                         beam=None, traj=False, end=1.8, dt=1.e-5, directory='./', Zlimit=-0.043, filename="Particles_Trajectories"):
   
    self.name = name
    self.date = date
    self.field = field
    self.source = source
    self.devices = devices
    self.beam = beam
    self.traj = traj
    self.dt = dt
    self.detector = detector
    self.Zlimit = Zlimit
    if not os.path.exists(directory):
      os.makedirs(directory)
    self.directory = directory
    self.filename = filename
    self.CalculateTrajectory(tStart=0, tEnd=end, dt=dt)
   
  def CalculateTrajectory(self, tStart, tEnd, dt ):
    """ Calculate the trajectories for all particles by integrating the equation of motion"""

    print("Calculating particles trajecories........")
    print("Running on "+str(cpu_count())+" cores")
    try:
      p = Pool()
      l = Manager().list()  # shared list to save the particles across the processes
      parallel_traj = partial(SingeParticleTrajectory, devices=self.devices, field=self.field, 
                              beam=self.beam, detector=self.detector, tStart=tStart, 
                              tEnd=tEnd, dt=dt, zlim=self.Zlimit, l=l, traj=self.traj)

      p.map(parallel_traj, self.source.particles)

    except Exception as e:
      print('Exception', e)

    finally:
      p.close()
      p.join()    
     
    self.SaveParticles(l, traj=self.traj)
    return True

  def SaveParticles(self, l, traj=False):
    """This function saves the initial and final particles' phase space 
       and temperatures to hdf5 files"""

    allParticles=[]
    for i in l:
      if i.reached:
        allParticles.append(i.FinalPhaseSpace)
    print(len(allParticles), "Particles Servived")
    f = h5py.File(self.directory+self.filename+".hdf5", "w")
    data = np.array(allParticles)
    f.create_dataset("particles", data=data)
    if traj:
      for p in l:
        f[f"trajectories/{p.ID}"] = np.array(p.trajectory)
      f.close()

  def LongestTrajectory(self):
    """ This function returns the longest trajectory of all particles 
        and used only for trajectories visualizations"""

    longest = 0
    for i in self.source.particles:
       if len(i.trajectory)>longest:
         longest = len(i.trajectory)
    return longest
