#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of CMInject
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this program for scientific work, you should correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

import argparse

import numpy as np
from cminject.definitions.detectors.simple_z import SimpleZDetector
from cminject.definitions.devices import FluidFlowFieldDevice
from cminject.definitions.devices.fluid_flow_field_device import FlowType
from cminject.definitions.particles.spherical import SphericalParticle
from cminject.definitions.property_updaters.brownian_motion import BrownianMotionPropertyUpdater
from cminject.definitions.setups import Setup
from cminject.definitions.sources.variable_distributions import VariableDistributionSource
from cminject.experiment import Experiment
from cminject.utils.args import SetupArgumentParser, dist_description


class FieldFreeTestSetup(Setup):
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        flow_file = args.flow_file
        dt = 1e-6

        devices = [
            FluidFlowFieldDevice(filename=flow_file, flow_type=FlowType.STOKES, offset=np.array([0, 0])),
            FluidFlowFieldDevice(filename=flow_file, flow_type=FlowType.STOKES, offset=np.array([0, 0.03]))
        ]

        detectors = [SimpleZDetector(i, 0.02 + 0.01*i) for i in range(3)]
        sources = [VariableDistributionSource(main_args.nof_particles, position=args.position, velocity=args.velocity,
                                              radius=args.radius, rho=args.density, subclass=SphericalParticle)]
        property_updaters = [BrownianMotionPropertyUpdater(field=devices[0].fields[0], dt=dt),
                             BrownianMotionPropertyUpdater(field=devices[1].fields[0], dt=dt)]

        # Construct the Experiment, still missing brownian motion
        exp = Experiment(
            detectors=detectors,
            sources=sources,
            devices=devices,
            property_updaters=property_updaters,
            time_step=1e-6,
            number_of_dimensions=2
        )

        return exp

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()

        parser.add_argument('-p', '--position',
                            help='Distribution description for the position.',
                            nargs='*', type=dist_description, required=True)
        parser.add_argument('-v', '--velocity',
                            help='Distribution description for the velocity.',
                            nargs='*', type=dist_description, required=True)
        parser.add_argument('-r', '--radius', help='Distribution description for the radius.',
                            type=dist_description, required=True)
        parser.add_argument('-rho', '--density', help='Density of the particles.',
                            type=float, required=True)
        parser.add_argument('-ff', '--flow-file', help='The flow field file (HDF5) for the first skimmer',
                            required=True, type=str)

        return parser

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
