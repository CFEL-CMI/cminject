from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad, dblquad
import numpy as np
from reader import *

class Fluid:
   """If the interaction field is EM field then the set_method function here will be used"""
   def __init__(self, density, dynamic_viscosity, temperature = 298, pressure=100, 
                  molar_mass=0.004002602 ,thermal_creep = 1, inflow_speed=20, outflow_pressure=0, method='LBM', directory='../../sandbox/',filename=None):
     self.density = density
     self.eta = dynamic_viscosity
     self.mu = dynamic_viscosity/density
     self.T = temperature
     self.k = thermal_creep
     self.P = pressure
     self.mM = molar_mass # for helium (default) 0.004002602 kilograms per mole
     self.inflow_speed =inflow_speed
     self.outflow_pressure = outflow_pressure
     self.method = method
     self.inletSpeed = inflow_speed
     self.directory = directory
     self.FlowField = False
     if method=='LBM' and filename is None:
       self.Run_LB_Code()
     if filename is not None:
       self.ReadFromFile(filename)

   def set_Boundary_conditions(self, inflow_v, outflow_p):
     self.inflow_v = inflow_v
     self.outflow_p = outflow_p
     
   def Run_LB_Code(self):
     from subprocess import call
     call(["/home/aminmuha/cmi-injector/bin/bgc",str(self.inletSpeed), str(self.eta), str(self.density), self.directory])


   def ReadFromFile(self, filename):
     x, y, z, Vx, Vy, Vz, P = ReadFromFile(filename)
     self.fvx = RegularGridInterpolator((x, y, z), Vx)
     self.fvy = RegularGridInterpolator((x, y, z), Vy)
     self.fvz = RegularGridInterpolator((x, y, z), Vz)
     self.fp = RegularGridInterpolator((x, y, z), P)
     self.FlowField = True
