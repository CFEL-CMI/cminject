class AerodynamicsLensStack(object):
    """ This class implement the aerodynamic lens stack. It is contains a several adl segments"""
    def __init__(self, position, segments):
        self.position = position
        self.segments = segments
    
    def add_fluid(self, fluid):
        self.fluid = fluid
    
    def set_boundary_condition(self, inflow, outflow):
        self.inflow = inflow
        self.outflow = outflow
        return

    def add_particle(self, particle):
        return particle

