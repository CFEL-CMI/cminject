import sys
sys.path.insert(0, '../lib')

from device.source import *
from device.als import *
from device.experiment import *
from field.interaction_field import *
from visualizer.visualizer3d import *
import matplotlib.pyplot as plt

# Create Particle
Radius = 1e-8 # 0.00000001
ParticleDensity = 100
NumberOfParticles = 6
SourceCoordinates = (0, 0, 0)
SigmaParticlesPosition = (0.003, 0.003, 0)
MuParticlesVelocity = (1, 0.5, -50)
SigmaParticlesVelocity = (0.5, 0.1, 0.5)

s = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=ParticleDensity  )

# Define aerodynamic lens stack
ADSPosition = (0, 0, 0)

# Each segment in the als is defined by the radius and the extention in the z-direction (R, Z)
Segments = [(0.01, -0.025), (0.001, -0.005), (0.01, -0.015), (0.001, -0.005)]

adl = AerodynamicsLensStack(ADSPosition, Segments)

# Define Fluid Object
FluidDensity = 0.016 
DynamicViscosity = 0.0012

f = Fluid( FluidDensity, DynamicViscosity, filename='3dfield.txt')


ExpName = 'Particles moving with constant velocity'
ExpDate = 'October'
SourceOfParticles = s
exp = Experiment( ExpName, ExpDate, SourceOfParticles, als=adl, field=f )

visualizer(exp)

