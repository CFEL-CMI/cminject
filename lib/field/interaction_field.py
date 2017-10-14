from scipy.interpolate import RegularGridInterpolator
import numpy as np

class interaction_field:
  """This class is represent the fields that are added to the experiment.
     The user may introduce electromagnetic of fluid field to the experiment"""
  
  def __init__(self, field_type=None, field=None):
    self.field_type = field_type
    self.field = field

class EM(interaction_field):
  """If the interaction field is EM field then the set_method function here will be used"""
  def __init__(self, intensity, energy):
     self.intensity = intensity
     self.energy = energy

  def set_method(self, method='imperical'):
     self.method = method
  
class Fluid(interaction_field):
   """If the interaction field is EM field then the set_method function here will be used"""
   def __init__(self, density, kinematic_viscosity, inflow_speed=0, outflow_pressure=0, method='LBM', filename=None):
     self.density = density
     self.mu = kinematic_viscosity
     self.inflow_speed =inflow_speed
     self.outflow_pressure = outflow_pressure
     self.method = method
     if filename is not None:
       self.ReadFromFile(filename)

   def set_Boundary_conditions(self, inflow_v, outflow_p):
     self.inflow_v = inflow_v
     self.outflow_p = outflow_p
     
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
