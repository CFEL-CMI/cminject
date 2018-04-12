from experiment.source import *
from experiment.device.als import *
from experiment.experiment import *
from field.interaction_field import *
from experiment.detector import Detector
import os

# Create Particle
Radius = 5e-7 # 0.00000001
ParticleDensity = 1050
NumberOfParticles = 3000
SourceCoordinates = (0, 0, -0.002)
SigmaParticlesPosition = (0.0005, 0.0005, 0)
MuParticlesVelocity = (0., 0., -50.)
SigmaParticlesVelocity = (0.0, 0.0, 0.0)

SourceOfParticles = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=ParticleDensity  )

# Define Fluid Object
FluidDensity = 0.012 
DynamicViscosity = 1.992e-5

f = Fluid( FluidDensity, DynamicViscosity, filename='BGC_50sccm.txt')  # read flow field from comsole file
d = Detector(-0.052)

directory='/Users/aminmuha/Documents/myproject/50sccm/0.0005mmParticlePositionInlet_3000Particles/'

ExpName = 'Buffer Gas Cell'
ExpDate = 'April2018'
exp = Experiment( ExpName, ExpDate, SourceOfParticles, end=0.8, field=f, detector=d, directory=directory )

