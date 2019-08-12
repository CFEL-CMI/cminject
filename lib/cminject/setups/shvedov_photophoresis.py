#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of CMInject.
# It is a reference implementation of a subset of possible experimental setups, together with command line argument
# parsing. Its purpose is to allow those simulations to be run without additional code needing to be written, and to
# showcase how the setup is created so people can extend it or write more complex setups.
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
from typing import Tuple, List

from cminject.definitions import Device, Field
from cminject.definitions.boundaries import CuboidBoundary
from cminject.definitions.detectors import SimpleZDetector
from cminject.definitions.fields.laser_field import ShvedovPhotophoreticLaserField
from cminject.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.definitions.sources import VariableDistributionSource
from cminject.experiment import Experiment
from cminject.setups.base import Setup
from cminject.utils.args import dist_description, SetupArgumentParser


class ShvedovVortexLaserDevice(Device):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

    def __init__(self):
        pp_field = ShvedovPhotophoreticLaserField(
            gas_temperature=293.15,
            gas_viscosity=1.76e-5,
            gas_thermal_conductivity=0.02546,  # W/(m*K) for nitrogen
            gas_density=1.161,  # mg/(cm^3) for nitrogen
            beam_power=0.01,  # W
            beam_waist_radius=8.4e-6,  # m
        )
        self.fields: List[Field] = [pp_field]
        self.boundary: CuboidBoundary = CuboidBoundary(
            intervals=[(-1e5, 1e5), (0, 0.02)]
        )
        super().__init__(fields=self.fields, boundary=self.boundary)


class ShvedovPhotophoresisSetup(Setup):
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        devices = [ShvedovVortexLaserDevice()]
        detectors = [SimpleZDetector(identifier=i, z_position=pos) for i, pos in enumerate(args.detectors)]
        sources = [VariableDistributionSource(main_args.nof_particles, position=args.position, velocity=args.velocity,
                                              radius=args.radius, rho=args.density,
                                              subclass=ThermallyConductiveSphericalParticle,
                                              thermal_conductivity=args.thermal_conductivity)]
        return Experiment(devices=devices, detectors=detectors, sources=sources, property_updaters=[],
                          time_interval=main_args.time_interval, delta_z_end=0.0, seed=main_args.seed,
                          number_of_dimensions=2)

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-d', '--detectors', help='The Z positions of the detectors', nargs='+',
                            type=float, required=True)
        parser.add_argument('-p', '--position',
                            help='Distribution description for the position.',
                            nargs='*', type=dist_description)
        parser.add_argument('-v', '--velocity',
                            help='Distribution description for the velocity.',
                            nargs='*', type=dist_description)
        parser.add_argument('-r', '--radius', help='Distribution description for the radius.',
                            type=dist_description)
        parser.add_argument('-rho', '--density', help='Density of the particles.',
                            type=float)
        parser.add_argument('-mu', '--thermal-conductivity', help='Thermal conductivity of the particles.',
                            type=float)

        parser.set_defaults(
            time_interval=[0.0, 1.0, 1e-5],
            thermal_conductivity=0.0266,
            rho=1775.0,
            density=1050.0,
            loglevel='warning',
            chunksize=1,  # same as None according to multiprocessing docs
        )
        return parser

    @staticmethod
    def validate_args(args: argparse.Namespace):
        # Verify dimensionality match for position description and dimensions parameter
        if len(args.position) != 2:
            raise argparse.ArgumentError("The length of the position description must be 2!")
        if len(args.velocity) != 2:
            raise argparse.ArgumentError("The length of the velocity description must be 2!")


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
