import sys
sys.path.insert(0, '../lib')

from device.source import *
from device.experiment import *
from visualizer.visualizer3d import *
import matplotlib.pyplot as plt


NumberOfParticles = 5
SourceCoordinates = (0, 0.01, 0)
SigmaParticlesPosition = (0, 0.003, 0)
MuParticlesVelocity = (10, 0, 0)
SigmaParticlesVelocity = (5, 0, 0)

s = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity  )

ExpName = 'Particles moving with constant velocity'
ExpDate = 'October'
SourceOfParticles = s
exp = Experiment( ExpName, ExpDate, SourceOfParticles)

visualizer(exp)

