import sys
sys.path.insert(0, '../lib')

from device.experiment import *
from device.asl import *
from device.detector import *
from device.source import *
from particle.particle import *
from field.interaction_field import *


# create a particle
radius = 0.00001     # m
temp = 293           # K
density = 10         # Kg/m3
thermal_conductivity = 1  
index_of_ref = 1.2 
particle = SphericalParticle( radius, temp, density, thermal_conductivity, index_of_ref )


devices=[]

# create a source
number_of_particles = 1000
sigma_of_particles_momentum = 3
source = Source(particle, number_of_particles, sigma_of_particles_momentum)



# create fluid object
FluidDensity = 0.0002      # Kg/m3
KinematicViscosity =0.12   # m2/s
InFlowSpeed = 10           # m/s
OutFlowPressure = 0        # Pascal
ModellingMethod = 'LBM'     # LBM/NS
helium = Fluid( FluidDensity, KinematicViscosity, InFlowSpeed, OutFlowPressure, ModellingMethod)


# create aerodynamic lens device
Hieght = 0.05  # m
RadiusInlet = 0.01 # m
RadiusLens = 0.005 # m
RadiusOutlet = 0.002 # m
adls = AerodynamicsLensDevice(Hieght, RadiusInlet, RadiusLens, RadiusOutlet, helium)

center = (0, 0, 100)
length = 0.5         # m
width =  0.5         # m
detector = Detector( center, length, width)

