from experiment.source import *
from experiment.device.als import *
from experiment.experiment import *
from field.interaction_field import *
from visualizer.visualizer3d import *
from experiment.detector import Detector
import os
import sys
import time
import datetime

sys.stdout = open(os.devnull, 'w') #supress print statements. use for multiprocessing

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M-%S')
filename = sys.argv[1] #flow-field file
vz = float(sys.argv[2]) #velocity in z direction. particle fly in -z
out = filename.split('.')[0] +'_'+st
outdir = sys.argv[3] #output directory
# Create Particle
Radius = 1.e-7 # 0.00000001 #in meters
ParticleDensity = 1050 #in kg/m^3
NumberOfParticles = int(sys.argv[4])
SourceCoordinates = (0, 0, 5.0000000000000E-4) #0,0,0 center front surface of inlet
SigmaParticlesPosition = (0.00025, 0.00025, 0.0000001) # normal distribution
MuParticlesVelocity = (0., 0., vz) #m/s
SigmaParticlesVelocity = (2.00, 2.00, 5.00) #

SourceOfParticles = Source( NumberOfParticles, SourceCoordinates , SigmaParticlesPosition, MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=ParticleDensity  )

ADSPosition = (0, 0, 0)

Segments = [(0.01, -0.025), (0.01, -0.005), (0.01, -0.015), (0.01, -0.01)]

adl = AerodynamicsLensStack(ADSPosition, Segments)


# Define Fluid Object
FluidDensity = 0.0009
DynamicViscosity = 1.02e-6

f = Fluid( FluidDensity, DynamicViscosity, filename=filename)  # read flow field from comsole file
d = Detector([-0.052])

filename = 'v'+str(vz).replace("-","")+'_'+out+'_'

ExpName = 'Buffer Gas Cell'
ExpDate = 'April2018'
exp = Experiment( ExpName, ExpDate, SourceOfParticles, end=1., dt=1.e-5, field=f, devices=[adl], detector=d, directory=outdir, filename=filename )
#visualizer(exp)
