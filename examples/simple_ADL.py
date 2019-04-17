# -*- coding: utf-8 -*-
#
# This file is part of CMInject


__doc__ = """Example input file demonstrating a simple aerodynamic-lens stack"""


from cminject.experiment.source import Source
from cminject.experiment.device.als import AerodynamicsLensStack
from cminject.experiment.experiment import Experiment
from cminject.experiment.detector import Detector

particlesource = Source(number_of_particles = 1,
                        position = (0,0,5e-4),
                        sigmaP = (0.0000001, 0.0000001, 0.0000001),
                        muV = (0., 4., -40.), 
                        sigmaV = (0., 0., 0.), radius=1.e-7, rho=1050.)


ADSPosition = (0, 0, 0)
Segments = [(0.01, -0.025), (0.001, -0.005), (0.01, -0.015), (0.001, -0.005)]
aerolensstack = AerodynamicsLensStack(ADSPosition, Segments)

devices = [aerolensstack,]

detector = Detector([-0.055]) 

epxeriment = Experiment('name', 'date', source=particlesource, devices=devices, detector=detector)



### Local Variables:
### mode: python
### fill-column: 80
### truncate-lines: t
### End:
