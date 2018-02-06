import sys
sys.path.insert(0, '../lib')

from device.source import *
from device.bgc import *
from device.experiment import *
from field.interaction_field import *
from visualizer.visualizer3d import *
import matplotlib.pyplot as plt

# Create Particle
Radius = 2e-7 # 0.00000001
ParticleDensity = 100
NumberOfParticles = 120
SourceCoordinates = (0.0005, 0.0155, 0.0155)
SigmaParticlesPosition = (0.0001, 0.0005, 0.0005)
MuParticlesVelocity = (15.0, 0., 0.)
SigmaParticlesVelocity = (0.5, 0., 0.)

s = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=ParticleDensity  )

# Define aerodynamic lens stack
BGCPosition = (0, 0.0155, 0.0155)


bgc = BufferGasCell(BGCPosition)

# Define Fluid Object
FluidDensity = 0.016 
DynamicViscosity = 0.00001992

f = Fluid( FluidDensity, DynamicViscosity, filename='bgcfield.txt')


ExpName = 'BufferGasCell'
ExpDate = 'Jan2018'
SourceOfParticles = s
exp = Experiment( ExpName, ExpDate, SourceOfParticles, devices=[bgc], field=f, end=30. )

visualizer(exp)

