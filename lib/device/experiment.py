from scipy.integrate import ode

class Experiment:
  """This class carry out the info about the experiment you want to simulate. 
     For example, the devices and the interaction fields you will include.
     The output files will have the name of and date of your experiment""" 

  def __init__(self, name, date, source, detector=None, adl=None, field=None):
    self.name = name
    self.date = date
    self.field = []
    self.device = []
    self.source = source
    self.CalculateTrajectory()

  def CalculateTrajectory(self, tStart=0, tEnd=2 ):
    """ Calculate the particle trajectories by integrating the equation of motion"""
    count = 0
    for i in self.source.particles:
      traj = 0
      integral = ode( i.get_v_and_a )
      integral.set_integrator('vode',method='BDF',with_jacobian=False,atol=1e-6,rtol=1e-6,first_step=1e-5,nsteps=1000)
      integral.set_initial_value( i.position + i.velocity, tStart )
      print "Calculating particle: ", count
      while integral.successful() and integral.t < tEnd:
        if traj ==1:
          i.trajectory.append( (integral.y[0], integral.y[1], integral.y[2]) )
          traj = 0
        integral.integrate(integral.t + 0.00001)
        traj+=1
      count+=1
    return integral
    

  def add_device(self, device):
    self.device.append(device)

  def add_field(self, device):
    self.field.append(field)

  def run(self):
   pass
