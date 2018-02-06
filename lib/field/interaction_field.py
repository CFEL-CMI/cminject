from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad, dblquad
import numpy as np

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
     if method=='LBM':
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
     f = open(filename)
     x=[]; y=[]; z=[]
     vx=[]; vy=[]; vz=[]
     p=[]
     for i in f:
       token = i.split()
       x.append(float(token[0]))
       y.append(float(token[1]))
       z.append(float(token[2]))
 
       if token[3]=='NaN':
         vx.append(0.0)
       else:
         vx.append(float(token[3]))

       if token[4]=='NaN':
         vy.append(0.0)
       else:
         vy.append(float(token[4]))

       if token[5]=='NaN':
         vz.append(0.0)
       else:
         vz.append(float(token[5]))

       if token[6]=='NaN':
         p.append(0.0)
       else:
         p.append(float(token[6]))

     x=sorted(set(x))
     y=sorted(set(y))
     z=sorted(set(z))
     n_x = len(x)
     n_y = len(y)
     n_z = len(z)
     c=0
     Vx = np.zeros((n_x, n_y, n_z))
     Vy = np.zeros((n_x, n_y, n_z))
     Vz = np.zeros((n_x, n_y, n_z))
     P = np.zeros((n_x, n_y, n_z))
     for i in range(n_z):
       for j in range(n_y):
         for k in range(n_x):
           Vx[k,j,i] = vx[c]
           Vy[k,j,i] = vy[c]
           Vz[k,j,i] = vz[c]
           P[k,j,i] = p[c]
           c+=1
     self.fvx = RegularGridInterpolator((x, y, z), Vx)
     self.fvy = RegularGridInterpolator((x, y, z), Vy)
     self.fvz = RegularGridInterpolator((x, y, z), Vz)
     self.fp = RegularGridInterpolator((x, y, z), P)
     self.FlowField = True
