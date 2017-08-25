from lib.device.experiment import *
from lib.device.asl import *
from lib.field.interaction_field import *

# start an experiment
exp = experiment('exp1', '28-08-2017')

# creat aerodynamic lens device
adls = als(0.05, 0.01, 0.005, 0.002, 1, 0.025)

# creat fluid object
helium = fluid(0.0002, 0.12)

# add boundary condition, inflow spead 5 m/s and pressure at the out let is 0
helium.set_Boundary_conditions(5, 0)

# add the fluid to the adls device
adls.add_fluid(helium, 'LBM')

# add device to the experiment
exp.add_device(adls) 


# run experiment
exp.run()
