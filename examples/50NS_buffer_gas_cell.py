from experiment.source import *
from experiment.device.als import *
from experiment.experiment import *
from field.interaction_field import *
from experiment.detector import Detector
import os

# Create Particle
Radius = 5e-7 # 0.00000001
ParticleDensity = 1050
NumberOfParticles = 1
SourceCoordinates = (0, 0, -0.000)
SigmaParticlesPosition = (0.000, 0.000, 0.000)
MuParticlesVelocity = (0., 0., -70.)
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
exp = Experiment( ExpName, ExpDate, SourceOfParticles, end=18., field=f, detector=d, directory=directory )

