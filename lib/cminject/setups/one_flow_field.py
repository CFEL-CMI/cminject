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
import logging
from typing import Tuple, List

from cminject.definitions import Device, Boundary
from cminject.definitions.boundaries import GridFieldBasedBoundary
from cminject.definitions.detectors import SimpleZDetector
from cminject.definitions.fields.fluid_flow_fields import StokesDragForceField, MolecularFlowDragForceField
from cminject.definitions.particles import SphericalParticle
from cminject.definitions.property_updaters import BrownianMotionPropertyUpdater, \
    BrownianMotionMolecularFlowPropertyUpdater
from cminject.definitions.sources import VariableDistributionSource
from cminject.experiment import Experiment
from cminject.setups.base import Setup
from cminject.utils.args import dist_description, SetupArgumentParser


class FluidFlowFieldDevice(Device):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

    def __init__(self, flow_type, *args, **kwargs):
        if flow_type == 'stokes':
            self.fields = [StokesDragForceField(*args, **kwargs)]
        elif flow_type == 'molecular_flow':
            self.fields = [MolecularFlowDragForceField(*args, **kwargs)]
        else:
            raise ValueError(
                f"Unknown flow type: {flow_type}. Available are: [stokes, molecular_flow]."
            )

        self.boundary: Boundary = GridFieldBasedBoundary(self.fields[0])
        super().__init__(fields=self.fields, boundary=self.boundary)


class OneFlowFieldSetup(Setup):
    """
    A setup for one fluid flow field, with a Stokes drag force and a molecular flow drag force available to choose
    from. Brownian motion can also be enabled for both models.
    """
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        t_start, t_end, dt = main_args.time_interval

        devices = [FluidFlowFieldDevice(flow_type=args.flow_type, filename=args.flow_field,
                                        temperature=args.flow_temperature,
                                        slip_correction_scale=args.flow_scale_slip or 1.0)]

        property_updaters = []
        if args.brownian:
            if args.flow_type == 'stokes':
                # noinspection PyTypeChecker
                property_updaters += [BrownianMotionPropertyUpdater(field=devices[0].fields[0], dt=dt)]
            elif args.flow_type == 'molecular_flow':
                # noinspection PyTypeChecker
                property_updaters += [BrownianMotionMolecularFlowPropertyUpdater(field=devices[0].fields[0], dt=dt)]

        detectors = [SimpleZDetector(identifier=i, z_position=pos) for i, pos in enumerate(args.detectors)]
        sources = [VariableDistributionSource(main_args.nof_particles, position=args.position, velocity=args.velocity,
                                              radius=args.radius, rho=args.density, subclass=SphericalParticle)]

        return Experiment(devices=devices, detectors=detectors, sources=sources, property_updaters=property_updaters,
                          time_interval=(t_start, t_end, dt), delta_z_end=0.001, seed=main_args.seed,
                          number_of_dimensions=args.dimensions)

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-f', '--flow-field', help='Flow field filename (hdf5 format)',
                            type=str, required=True)
        parser.add_argument('-F', '--flow-type', help='The type of flow model to use.',
                            type=str, choices=['stokes', 'molecular_flow'], default='stokes')
        parser.add_argument('-D', '--dimensions', help='# of spatial dimensions',
                            type=int)
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

        parser.add_argument('-d', '--detectors', help='The Z positions of the detectors', nargs='+',
                            type=float, required=True)
        parser.add_argument('-B', '--brownian', help='Enable Brownian motion', action='store_true')

        parser.add_argument('-ft', '--flow-temperature',
                            help='Temperature of the gas in the flow field. NOTE: This should be close or equal to the '
                                 'temperature the flow field file was created with, unless you are sure that the used '
                                 'model and assumed parameters work across a large range of temperatures',
                            type=float)
        parser.add_argument('-fs', '--flow-scale-slip', help='Slip scale factor of the flow field',
                            type=float)

        parser.set_defaults(
            dimensions=3,
            radius={'kind': 'gaussian', 'mu': 2.45e-7, 'sigma': 7.5e-9},
            position=[
                {'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.0002},
                {'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.0002},
                {'kind': 'gaussian', 'mu': 0.0048, 'sigma': 0.00001}
            ],
            velocity=[
                {'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.1},
                {'kind': 'gaussian', 'mu': 0.7, 'sigma': 0.3},
                {'kind': 'gaussian', 'mu': -43.0, 'sigma': 2.0},
            ],
            density=1050.0,
        )

        return parser

    @staticmethod
    def validate_args(args: argparse.Namespace):
        # Verify dimensionality match for position description and dimensions parameter
        if len(args.position) != args.dimensions:
            raise argparse.ArgumentError(
                "The length of the position description vectors must match the simulation dimensionality!"
            )
        if len(args.velocity) != args.dimensions:
            raise argparse.ArgumentError(
                "The length of the velocity description vectors must match the simulation dimensionality!"
            )


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
