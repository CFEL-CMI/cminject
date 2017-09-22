import sys
sys.path.insert(0, '../lib')

from device.experiment import *
from device.asl import *
from device.detector import *
from device.source import *
from particle.particle import *
from field.interaction_field import *


# create a particle
Radius = 0.00001     # m
Temp = 293           # K
Density = 10         # Kg/m3
ThermalConductivity = 1  
IndexOfRef = 1.2 
particle = SphericalParticle( Radius, Temp, Density, ThermalConductivity, IndexOfRef )


devices=[]

# create a source
NumberOfParticles = 1000
SigmaOfParticlesMomentum = 3
source = Source(particle, NumberOfParticles, SigmaOfParticlesMomentum)



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

# create a laser field
Intensity = 100000  # W/cm2
Energy =  300       # eV
field = EM( Intensity, Energy)

# create a detector
center = (0, 0, 100)
length = 0.5         # m
width =  0.5         # m
detector = Detector( center, length, width)

# create an experiment, your output files will be named after the experiment name and date
exp = Experiment('Exp1', '22-9-2017', detector, source, adl=[adls], field=[field]) 
