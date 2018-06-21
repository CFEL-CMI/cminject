from experiment.source import *
from experiment.device.als import *
from experiment.experiment import *
from field.interaction_field import *
from experiment.detector import Detector
#from visualizer.visualizer3d import *
from experiment.device.als import *
import os
import sys

filename = sys.argv[1]
vz = float(sys.argv[2])
out = filename.split('.')[0]
outdir = sys.argv[3]
# Create Particle
Radius = 1.e-7 # 0.00000001
ParticleDensity = 1050
NumberOfParticles = 100
SourceCoordinates = (0, 0, 0.0065)
SigmaParticlesPosition = (0.0000001, 0.0000001, 0.0000001)
MuParticlesVelocity = (0., 0., vz)
SigmaParticlesVelocity = (.1, .1, 0.1)

SourceOfParticles = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=ParticleDensity  )

ADSPosition = (0, 0, 0)

Segments = [(0.01, -0.025), (0.01, -0.005), (0.01, -0.015), (0.01, -0.0055)]

adl = AerodynamicsLensStack(ADSPosition, Segments)


# Define Fluid Object
FluidDensity = 0.0009 
DynamicViscosity = 1.02e-6

f = Fluid( FluidDensity, DynamicViscosity, filename=filename)  # read flow field from comsole file
d = Detector([-0.01, -0.052])

filename = 'v'+str(vz).replace("-","")+'_'+out+'_'

ExpName = 'Buffer Gas Cell'
ExpDate = 'April2018'
exp = Experiment( ExpName, ExpDate, SourceOfParticles, end=1., dt=1.e-5, field=f, devices=[adl], detector=d, directory=outdir, filename=filename )
#visualizer(exp)
