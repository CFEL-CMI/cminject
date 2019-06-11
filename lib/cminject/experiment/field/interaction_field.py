# -*- coding: utf-8 -*-
import sys
import os
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad, dblquad
import numpy as np
from math import exp, sqrt, pi
from cminject.experiment.field.reader import *
from cminject.experiment_refactored.tools.structured_txt_hdf5_tools import hdf5_to_data_grid

class Fluid:
   """If the interaction field is EM field then the set_method function here will be used"""
   def __init__(self, density, dynamic_viscosity, temperature = 4.0, 
                  pressure=100, molar_mass=0.004002602 ,thermal_creep = 1, inflow_speed=20, 
                    outflow_pressure=0, Kn=912.0, kinetic_d = 260.e-12, specific_gas_const = 2077., 
                      molar_heat_capacity = 12.5, method='LBM', conv=0.00001, lattice_velocity=0.01, 
                         specific_heat_capacity = 3116., mGas=6.6e-27, mGasMol=0.004002602,  directory='./',
                            speed_of_sound = 117.7, ScaleSlip=1., filename=None, new=True):

     self.density = density
     self.eta = dynamic_viscosity/density #kinematic_viscosity
     self.mu = dynamic_viscosity  #kinematic_viscosity * density
     self.T = temperature
     self.k = thermal_creep
     self.P = pressure
     self.mM = molar_mass # for helium (default) 0.004002602 kilograms per mole
     self.inflow_speed =inflow_speed
     self.outflow_pressure = outflow_pressure
     self.Kn = Kn  
     self.cv = molar_heat_capacity
     self.kinetic_d = kinetic_d
     self.mGas = mGas
     self.Rs = specific_gas_const
     self.cp = specific_heat_capacity
     self.mGasMol = mGasMol
     self.c = speed_of_sound
     self.SlipCorrection = 1 + Kn * (1.2310 + (0.4695 * exp(-1.1783/Kn)))
     self.method = method
     self.inletSpeed = inflow_speed
     self.directory = directory
     self.conv = conv
     self.U = lattice_velocity
     self.ScaleSlip = ScaleSlip
     self.new = new
     self.FlowField = False
     self.pressure = 0
     self.vel = (0,0,0)

     if method=='LBM' and filename is None:
       print("Running Lattice Boltzmann Code")
       self.Run_LB_Code()
     if filename is not None:
       self.ReadFromFile(filename)

   def set_Boundary_conditions(self, inflow_v, outflow_p):
     self.inflow_v = inflow_v
     self.outflow_p = outflow_p
     
   def Run_LB_Code(self):
     from subprocess import call
     if self.new:
       call([PATH+"bgc",str(self.inletSpeed), str(self.eta), 
               str(self.density), str(self.conv), str(self.U), self.directory])
     self.ReadFromFile(self.directory+"/LBM/vtkData/data/", True)     

   def MeanFreePath(self):
     """ Return mean free path of the gas"""
     k = 1.38065e-23     
     lamda = (self.mu/self.pressure)*sqrt(pi*k*self.T/(2*self.mGas))
     return lamda

   def ThermalConductivity(self):
     """ returns thermal conductivity based on ideal gas law
         http://hyperphysics.phy-astr.gsu.edu/hbase/thermo/thercond.html"""

     v = sqrt(self.vel[0]**2 + self.vel[1]**2 + self.vel[2]**2)
     l = self.MeanFreePath()
     k = self.pressure * v * l * self.cv/(3 * self.mM * self.Rs * self.T)
#     print "pressure", self.pressure, "vel", v, "MFP", l, "thermal C", k
     return k

   def Pr(self):
     """This function calculates Prandtl number"""
     return self.cp * self.mu / self.ThermalConductivity()

   def Mach(self):
     """This function calculates Mach number"""
     v = sqrt(self.vel[0]**2 + self.vel[1]**2 + self.vel[2]**2)
     return v/self.c
 
   def Reynolds(self, D):
     """ Calculates Reynold number given characteristic length D"""
     rho = self.pressure /(self.T * self.Rs)
     v = sqrt(self.vel[0]**2 + self.vel[1]**2 + self.vel[2]**2)
     return rho * v * D / self.mu

   def Nusselt(self, D, muS):
     """ Calculates Nusselt number given muS, which is the dynamic 
         viscosity at the surface temperature of nanoparticle"""
     Re = self.Reynolds(D)
     Pr = self.Pr()
     Nu = 2 + (0.4*Re**0.5 + 0.06*Re**0.67) * (Pr**0.4) * (self.mu/muS)**0.25
     return Nu

   def h(self, Nu, D):
     """The convective heat transfer coefficient h strongly depends on
         the fluid properties and roughness of the solid surface, and 
         the type of the fluid flow (laminar or turbulent)"""
     h = Nu * self.ThermalConductivity() / D
     return h
 
   def ReadFromFile(self, filename):
     xyz,data_grid = hdf5_to_data_grid(filename)
     x,y,z = xyz
     self.fdrag = RegularGridInterpolator((x, y, z), data_grid)
     self.FlowField = True
