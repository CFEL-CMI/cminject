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

from cminject.boundaries.grid_field_based import GridFieldBasedBoundary
from cminject.detectors.simple_z import SimpleZDetector
from cminject.devices.fluid_flow_device import FluidFlowDevice, FlowType
from cminject.fields.fluid_flow import StokesDragForceField, MolecularFlowDragForceField
from cminject.particles.spherical import SphericalParticle, ThermallyConductiveSphericalParticle
from cminject.actions.brownian_motion import\
    StokesBrownianMotionStep, MolecularFlowBrownianMotionStep
from cminject.base import Setup, Device
from cminject.sources.variable_distributions import VariableDistributionSource

from cminject.experiment import Experiment
from cminject.utils.args import dist_description, SetupArgumentParser


class OneFlowFieldSetup(Setup):
    """
    A setup for one fluid flow field, with a Stokes drag force and a molecular flow drag force available to choose
    from. Brownian motion can also be enabled for both models.
    """
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        t_start, t_end = main_args.time_interval
        dt = main_args.time_step
        source_kwargs = {
            'number_of_particles': main_args.nof_particles,
            'position': args.position, 'velocity': args.velocity,
            'radius': args.radius, 'density': args.density
        }
        ff_kwargs = {
            'filename': args.flow_field, 'brownian_motion': args.brownian,
            'flow_type': FlowType.MOLECULAR_FLOW if args.flow_type == 'molecular_flow' else FlowType.STOKES,
            'temperature': args.flow_temperature
        }
        if args.flow_type == 'molecular_flow':
            # Need to set subclass and pass starting temperature when using molecular flow
            source_kwargs = {
                **source_kwargs, 'subclass': ThermallyConductiveSphericalParticle,
                'subclass_kwargs': {'temperature': args.particle_temperature}
            }
        else:
            # Need to pass along slip correction scale if using stokes flow
            ff_kwargs = {**ff_kwargs, 'slip_correction_scale': args.flow_scale_slip or 1.0}

        # Construct Experiment
        experiment = Experiment(number_of_dimensions=args.dimensions, time_interval=(t_start, t_end), time_step=dt,
                                random_seed=main_args.seed)
        # Add Source and Device
        experiment.add_source(VariableDistributionSource(**source_kwargs))
        experiment.add_device(FluidFlowDevice(**ff_kwargs))
        # Simple Z detectors at every given position
        for i, pos in enumerate(args.detectors or []):
            experiment.add_detector(SimpleZDetector(identifier=i, z_position=pos))
        return experiment

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-f', '--flow-field', help='Flow field filename (hdf5 format)',
                            type=str, required=True)
        parser.add_argument('-F', '--flow-type', help='The type of flow model to use.',
                            type=str, choices=['stokes', 'molecular_flow'], default='stokes')
        parser.add_argument('-D', '--dimensions', help='# of spatial dimensions',
                            type=int, required=True)
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

        parser.add_argument('-d', '--detectors', help='The Z positions of the detectors', nargs='*',
                            type=float, required=False)
        parser.add_argument('-B', '--brownian', help='Enable Brownian motion', action='store_true')

        parser.add_argument('-ft', '--flow-temperature',
                            help='Temperature of the gas in the flow field. NOTE: This should be close or equal to the '
                                 'temperature the flow field file was created with, unless you are sure that the used '
                                 'model and assumed parameters work across a large range of temperatures',
                            type=float)
        parser.add_argument('-fs', '--flow-scale-slip', help='Slip scale factor of the flow field',
                            type=float)

        parser.add_argument('-pt', '--particle-temperature',
                            help='(Initial) temperature of the particles in K. Ignored if --flow-type=stokes, '
                                 'required if --flow-type=molecular_flow.',
                            type=float, required=False)
        return parser

    @staticmethod
    def validate_args(args: argparse.Namespace):
        # Verify dimensionality match for position description and dimensions parameter
        if len(args.position) != args.dimensions:
            raise ValueError(
                "The length of the position description vectors must match the simulation dimensionality!"
            )
        if len(args.velocity) != args.dimensions:
            raise ValueError(
                "The length of the velocity description vectors must match the simulation dimensionality!"
            )
        # Verify that particle temperature is passed when using molecular flow
        if args.flow_type == 'molecular_flow' and args.particle_temperature is None:
            raise ValueError(
                "--particle-temperature must be passed when --flow-type=molecular_flow!"
            )


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
