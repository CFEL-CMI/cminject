from experiment.source import *
from experiment.device.bgc import *
from experiment.experiment import *
from field.interaction_field import *
import time
# Create Particle
Radius = 2e-7 # 0.00000001
ParticleDensity = 100
NumberOfParticles = 10
SourceCoordinates = (0.005, 0.0154, 0.0154)
SigmaParticlesPosition = (0.0001, 0.00001, 0.0001)
MuParticlesVelocity = (25.0, 2., 0.5)
SigmaParticlesVelocity = (1., 0., 0.)

s = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=ParticleDensity  )

# Define aerodynamic lens stack
BGCPosition = (0, 0.0154, 0.0154)


bgc = BufferGasCell(BGCPosition)

# Define Fluid Object
FluidDensity = 0.00066
KinematicViscosity = 0.03
directory = "/Users/aminmuha/Documents/myproject/30ms/"
f = Fluid( FluidDensity, KinematicViscosity, method='LBM', new=False, inflow_speed=20., directory=directory)

ExpName = 'BufferGasCell'
ExpDate = '10Feb2018'
SourceOfParticles = s
t1=time.time()
exp = Experiment( ExpName, ExpDate, SourceOfParticles, devices=[bgc], field=f, end=0.1, directory=directory, dt=(1.e-5), traj=1 )
#visualizer(exp)
print time.time()-t1
