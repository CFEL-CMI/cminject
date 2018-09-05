# -*- coding: utf-8 -*-
#
# This file is part of CMInject


import cminject.apparatus.apparatus



class als(apparatus.apparatus):
    """ This class implement the aerodynamic lens stack. It is contains a several adl segments"""
    def __init__(self, no_of_segments, seperation, adl):
        super(als, self).__init__()
        self.no_of_segments = no_of_segments
        self.adl = adl
        self.seperation = seperation


    def add_fluid(self, fluid):
        self.fluid = fluid


    def set_boundary_condition(self, inflow, outflow):
        self.inflow = inflow
        self.outflow = outflow
        return


    def add_particle(self, particle):
        return particle



class segment:
    """This class implement the aerodynamic lens. The constructor take the width and the length of each segment"""
    def __init__(self, length, width):
        self.length = length
        self.width = width



### Local Variables:
### fill-column: 80
### truncate-lines: t
### End:
