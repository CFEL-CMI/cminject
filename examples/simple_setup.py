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

"""
This file defines a simple example setup. It is, at its core, a reduced version of
cminject.definitions.setups.OneFlowFieldSetup, and should only be used as an example to derive own setups from.

See README.md for further information on how to run this setup.
"""

# We need argparse for parsing arguments from the command line
import argparse

# Import the most basic definitions for defining a Setup which returns an Experiment
# and has arguments parsed by a SetupArgumentParser
from cminject.experiment import Experiment
from cminject.definitions.setups import Setup
from cminject.utils.args import SetupArgumentParser

# Import some concrete class implementations to be able to define a simple setup:
# A source which can generate different distributions
from cminject.definitions.sources import VariableDistributionSource
# A detector which is positioned at some Z position
from cminject.definitions.detectors import SimpleZDetector
# A device for fluid flow based on an HDF5 file, together with an enumeration of the possible models to pick (FlowType).
# (this setup only imports FlowType to access FlowType.STOKES, and doesn't allow picking by a user)
from cminject.definitions.devices.fluid_flow_field_device import FluidFlowFieldDevice, FlowType
# A PropertyUpdater for modeling Brownian motion based on a Stokes drag force.
from cminject.definitions.property_updaters import BrownianMotionPropertyUpdater


class SimpleSetup(Setup):
    """
    A simple setup which will simulate a single flow field with particles starting from a fixed position and velocity
    distribution and radius. The only parameters available are the flow field filename and the density (rho) of the
    particle material. The density is there just to show how additional arguments can be added.
    """
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        device = FluidFlowFieldDevice(filename=args.filename, flow_type=FlowType.STOKES)
        dt = 1e-5

        return Experiment(
            detectors=[
                SimpleZDetector(identifier=0, z_position=0.0),
                SimpleZDetector(identifier=1, z_position=0.01),
                SimpleZDetector(identifier=2, z_position=0.02),
            ],
            sources=[VariableDistributionSource(
                number_of_particles=main_args.nof_particles,
                position=[
                    {'kind': 'radial_gaussian', 'mu': 0.0, 'sigma': 1e-3},
                    0.0
                ],
                velocity=[
                    {'kind': 'gaussian', 'mu': 0.0, 'sigma': 1e-3},
                    {'kind': 'gaussian', 'mu': 2.0, 'sigma': 0.5}
                ],
                radius=5e-9,
                rho=args.rho
            )],
            devices=[device],
            property_updaters=[
                BrownianMotionPropertyUpdater(field=device.fields[0], dt=dt)
            ],
            time_step=dt,
            number_of_dimensions=2,
        )

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-f', '--filename', help='The filename of the flow field (HDF5).', type=str)
        parser.add_argument('--rho', help='The density of the particle material [kg/m^3].', type=float, default=1050.0)
        return parser

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
