import sys
sys.path.insert(0, '../lib')

from device.source import *
from device.experiment import *
from visualizer.visualizer import *
import matplotlib.pyplot as plt


NumberOfParticles = 10
SourceCoordinates = (0, 0.01)
SigmaParticlesPosition = 0.003
MuParticlesVelocity = 10
SigmaParticlesVelocity = 5

s = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity  )

ExpName = 'Particles moving with constant velocity'
ExpDate = 'October'
SourceOfParticles = s
exp = Experiment( ExpName, ExpDate, SourceOfParticles)

visualizer(exp)

