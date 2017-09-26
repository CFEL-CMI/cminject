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
   def __init__(self, density, kinematic_viscosity, inflow_speed, outflow_pressure, method='LBM'):
     self.density = density
     self.kinematic_viscosity = kinematic_viscosity
     self.inflow_speed =inflow_speed
     self.outflow_pressure = outflow_pressure
     self.method = method

   def set_method(self, method='LB'):
     self.method = method

   def set_Boundary_conditions(self, inflow_v, outflow_p):
     self.inflow_v = inflow_v
     self.outflow_p = outflow_p

   def ReadFromFile(self, filename):
     rawdata = open(filename,'r')
     # Read the file contents and generate a list with each line
     lines = rawdata.readlines()
     rawdata.close()
     #define variables
     rField=[]
     zField=[]
     vrField=[]
     vzField=[]
     pField=[]
     #write the lines into variables
     for line in lines:
        #remove eol
        newline=line[:-1]
        linesplit=newline.split(";")
        rField.append(float(linesplit[0]))
        zField.append(float(linesplit[1]))
        vrField.append(float(linesplit[2]))
        vzField.append(float(linesplit[3]))
        pField.append(float(linesplit[4]))
     rField=np.array(rField)
     zField=np.array(zField)
     vrField=np.array(vrField)
     vzField=np.array(vzField)
     pField=np.array(pField)
     return rField, zField, vrField, vzField, pField
    
