class AerodynamicsLensStack(object):
    """ This class implement the aerodynamic lens stack. It is contains a several adl segments"""
    def __init__(self, position, segments):
        self.position = position
        self.segments = segments
    
    def ParticleInside(self, position):
          x, y, z = position
          inside = False
          offset = self.position[2]
          for j in self.segments:
             if (x**2 +y**2) < j[0]**2 and abs(z) >= abs(offset) and abs(z) < abs(offset+j[1]):
               return True
             offset = offset + j[1]
          return False

    def add_fluid(self, fluid):
        self.fluid = fluid
    
    def set_boundary_condition(self, inflow, outflow):
        self.inflow = inflow
        self.outflow = outflow
        return

    def add_particle(self, particle):
        return particle

