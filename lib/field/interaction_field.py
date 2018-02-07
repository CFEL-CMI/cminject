import sys
import os
PATH = os.environ["CMIPATH"]+"bin/"
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad, dblquad
import numpy as np
from reader import *

class Fluid:
   """If the interaction field is EM field then the set_method function here will be used"""
   def __init__(self, density, dynamic_viscosity, temperature = 298, pressure=100, 
                  molar_mass=0.004002602 ,thermal_creep = 1, inflow_speed=20, outflow_pressure=0, 
                  method='LBM', conv=0.01, directory=os.environ["CMIPATH"]+'/sandbox/',filename=None):

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
     self.conv = conv
     self.FlowField = False
     if method=='LBM' and filename is None:
       print "Running Lattice Boltzmann Code"
       self.Run_LB_Code()
     if filename is not None:
       self.ReadFromFile(filename)

   def set_Boundary_conditions(self, inflow_v, outflow_p):
     self.inflow_v = inflow_v
     self.outflow_p = outflow_p
     
   def Run_LB_Code(self):
     from subprocess import call
     call([PATH+"bgc",str(self.inletSpeed), str(self.eta), str(self.density), str(self.conv), self.directory])
     self.ReadFromFile(self.directory+"/LBM/vtkData/data/", True)     

   def ReadFromFile(self, filename, vtk=False):
     if vtk:
       x, y, z, Vx, Vy, Vz, P = ReadVTK(filename)
     else:
       x, y, z, Vx, Vy, Vz, P = ReadText(filename)
     self.fvx = RegularGridInterpolator((x, y, z), Vx)
     self.fvy = RegularGridInterpolator((x, y, z), Vy)
     self.fvz = RegularGridInterpolator((x, y, z), Vz)
     self.fp = RegularGridInterpolator((x, y, z), P)
     self.FlowField = True
