# -*- coding: utf-8 -*-
#
# This file is part of CMInject


__doc__ = """Example input file demonstrating a simple aerodynamic-lens stack"""


from cminject.experiment.source import *
from cminject.experiment.device.als import *
from cminject.experiment.experiment import *


particlesource = Source(number_of_particles = int(1e6),
                        position = (0,0,5e-4),
                        sigmaP = (0.0000001, 0.0000001, 0.0000001),
                        MuParticlesVelocity, SigmaParticlesVelocity, radius=Radius, rho=ParticleDensity)

fluid = Fluid(flowfield = 'filename',
              fluid_density = 0.0009,
              dynamic_viscosity = 1.02e-6)


aerolensstack = ALS(origin=(0,0,0),
                    fluid = fluid)


devices = [aerolensstack,]

detector = ...

epxeriment = Experiment(particlesource, devices, detector)

#might experiment.run()






### Local Variables:
### mode: python
### fill-column: 80
### truncate-lines: t
### End:
