import sys
sys.path.insert(0, '../lib')

from device.source import *
from device.experiment import *
from field.interaction_field import *
from visualizer.visualizer3d import *
import matplotlib.pyplot as plt

# Create Particle
radius = 0.00001
density = 100
NumberOfParticles = 5
SourceCoordinates = (0, 0, 0)
SigmaParticlesPosition = (0.002, 0.003, 0)
MuParticlesVelocity = (0, 0, -5)
SigmaParticlesVelocity = (0, 0, 2)

s = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=radius, rho=density  )

f = Fluid( 0.01, 0.01, 10, 0, filename='3dfield.txt')


ExpName = 'Particles moving with constant velocity'
ExpDate = 'October'
SourceOfParticles = s
exp = Experiment( ExpName, ExpDate, SourceOfParticles, field=f )

visualizer(exp)

