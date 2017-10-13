import sys
sys.path.insert(0, '../lib')

from device.source import *
from device.als import *
from device.experiment import *
from field.interaction_field import *
from visualizer.visualizer3d import *
import matplotlib.pyplot as plt

# Create Particle
Radius = 0.00001
density = 100
NumberOfParticles = 4
SourceCoordinates = (0, 0, 0)
SigmaParticlesPosition = (0.003, 0.003, 0)
MuParticlesVelocity = (1, 0.5, -5)
SigmaParticlesVelocity = (0.5, 0.1, 0.5)

s = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=density  )

# Define aerodynamic lens stack
ADSPosition = (0, 0, 0)
Segments = [(0.01, -0.025), (0.001, -0.005), (0.01, -0.015), (0.001, -0.005)]
adl = AerodynamicsLensStack(ADSPosition, Segments)


f = Fluid( 0.01, 0.01, 10, 0, filename='3dfield.txt')


ExpName = 'Particles moving with constant velocity'
ExpDate = 'October'
SourceOfParticles = s
exp = Experiment( ExpName, ExpDate, SourceOfParticles, als=adl, field=f )

visualizer(exp)

