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
  


class fluid(interaction_field):
  """If the interaction field is EM field then the set_method function here will be used"""
   def __init__(self, density, kinematic_viscosity, inflow_velocity):
     self.density = density
     self.kinematic_viscosity = kinematic_viscosity
     self.inflow_velocity = inflow_velocity

   def set_method(self, method='LB'):
     self.method = method

   def set_Boundary_conditions(self, inflow_v, outflow_p):
     self.inflow_v = inflow_v
     self.outflow_p = outflow_p


    
