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

See README.md in this directory for further information on how to run this setup.
"""

# We need argparse for parsing arguments from the command line
import argparse

# Import the most basic definitions for defining a Setup which returns an Experiment
# and has arguments parsed by a SetupArgumentParser
from cminject.experiment import Experiment
from cminject.base import Setup
from cminject.utils.args import SetupArgumentParser
# Import some concrete class implementations to be able to define a simple setup:
# * A source which can generate different distributions
from cminject.sources import VariableDistributionSource
# * A detector which is positioned at some Z position
from cminject.detectors import SimpleZDetector
# * A device for fluid flow based on an HDF5 file, together with FlowType, an enumeration of the possible models to use
#   (this setup only imports FlowType to access FlowType.STOKES, and doesn't allow picking by a user)
from cminject.devices.fluid_flow_device import FluidFlowDevice, FlowType


class SimpleSetup(Setup):
    """
    A simple setup which will simulate a single flow field with particles starting from a fixed position and velocity
    distribution and radius. The only parameters available are the flow field filename and the density (rho) of the
    particle material. The density is there just to show how additional arguments can be added.
    """
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        dt = 1e-5
        experiment = Experiment(number_of_dimensions=2, time_step=dt, time_interval=(0, 1))
        experiment.add_source(VariableDistributionSource(
            number_of_particles=main_args.nof_particles,
            position=[{'kind': 'radial_gaussian', 'mu': 0.0, 'sigma': 1e-3}, 0.0],
            velocity=[{'kind': 'gaussian', 'mu': 0.0, 'sigma': 1e-3}, {'kind': 'gaussian', 'mu': 2.0, 'sigma': 0.5}],
            radius=5e-9,
            density=args.rho
        ))
        experiment.add_device(FluidFlowDevice(filename=args.filename, flow_type=FlowType.STOKES, brownian_motion=True))
        for i, z in enumerate([0.0, 0.01, 0.02]):
            experiment.add_detector(SimpleZDetector(z, z))

        return experiment

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-f', '--filename', help='The filename of the flow field (HDF5).', type=str, required=True)
        parser.add_argument('--rho', help='The density of the particle material [kg/m^3].', type=float, default=1050.0)
        return parser


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
