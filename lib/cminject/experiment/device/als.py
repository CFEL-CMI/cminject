class AerodynamicsLensStack(object):
    """ This class implement the aerodynamic lens stack. It is contains a several adl segments"""
    def __init__(self, position, segments):
        self.position = position
        self.segments = segments
        self.minz = position[2]
        self.maxz = position[2]
        for i in segments: 
          self.maxz += abs(i[1])
 
    def ParticleInside(self, position):
          """ This function take the position of the particle and return true if the particle inside the als"""
          x, y, z = position
          inside = False
          offset = self.position[2]
          for j in self.segments:
             if (x**2 +y**2) < j[0]**2 and abs(z) >= abs(offset) and abs(z) < abs(offset+j[1]):
               return True
             offset = offset + j[1]
          #print("The particle is outside the ADL, the final coordinates are", position)
          return False

    def add_fluid(self, fluid):
        self.fluid = fluid
    
    def set_boundary_condition(self, inflow, outflow):
        self.inflow = inflow
        self.outflow = outflow
        return

    def add_particle(self, particle):
        return particle

