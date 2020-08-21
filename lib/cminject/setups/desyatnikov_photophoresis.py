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
Setups for simulating photophoretic effects with the model described in: A. Desyatnikov, V. Shvedov, A. Rode,
W. Krolikowski, and Y. Kivshar, "Photophoretic manipulation of absorbing aerosol particles with vortex beams: theory
versus experiment," Opt. Express  17, 8201-8211 (2009).
"""

import argparse

import numpy as np

from cminject.base import Setup, Device
from cminject.boundaries.simple import InfiniteBoundary
from cminject.detectors.simple_z import SimpleZDetector
from cminject.devices.desyatnikov_photophoresis_device import DesyatnikovPhotophoresisDevice
from cminject.experiment import Experiment
from cminject.fields.function_based import FunctionField
from cminject.particles.spherical import ThermallyConductiveSphericalParticle
from cminject.sources.variable_distributions import VariableDistributionSource
from cminject.utils.args import dist_description, SetupArgumentParser


class DesyatnikovPhotophoresisSetup(Setup):
    """
    A setup simulating photophoretic effects with the model described in: A. Desyatnikov, V. Shvedov, A. Rode,
    W. Krolikowski, and Y. Kivshar, "Photophoretic manipulation of absorbing aerosol particles with vortex beams: theory
    versus experiment," Opt. Express  17, 8201-8211 (2009).
    """
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        experiment = Experiment(time_interval=main_args.time_interval, time_step=main_args.time_step,
                                random_seed=main_args.seed, number_of_dimensions=2)
        experiment.add_source(VariableDistributionSource(
            main_args.nof_particles, position=args.position, velocity=args.velocity, radius=args.radius,
            density=args.rho, subclass=ThermallyConductiveSphericalParticle,
            subclass_kwargs={'specific_heat': 0.0, 'temperature': 293.15,
                             'thermal_conductivity': args.gas_thermal_conductivity}
        ))
        experiment.add_device(DesyatnikovPhotophoresisDevice(
            r_boundary=args.boundary[:2], z_boundary=args.boundary[2:],
            gas_temperature=args.gas_temperature, gas_viscosity=args.gas_viscosity, gas_mass=args.gas_mass,
            gas_thermal_conductivity=args.gas_thermal_conductivity, gas_density=args.gas_density,
            beam_power=args.beam_power, beam_waist_radius=args.beam_waist_radius
        ))
        if args.gravity:
            experiment.add_device(
                Device(boundary=InfiniteBoundary(), fields=[FunctionField(lambda p, t: np.array([0, -9.81]))])
            )
        for i, pos in enumerate(args.detectors):
            experiment.add_detector(SimpleZDetector(identifier=i, z_position=pos))

        experiment.z_boundary = tuple(args.boundary[2:])
        return experiment

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()
        add_arg = lambda *args, **kwargs: parser.add_argument(*args, **kwargs)  # Shorthand functions for readability
        add_float_arg = lambda *args, **kwargs: add_arg(*args, **kwargs, type=float)
        add_dist_arg = lambda *args, **kwargs: add_arg(*args, **kwargs, type=dist_description)

        add_float_arg('-d', '--detectors', help='The Z positions of the detectors', nargs='+', required=True)
        add_dist_arg('-p', '--position', help='Distribution description for the position.', nargs='*')
        add_dist_arg('-v', '--velocity', help='Distribution description for the velocity.', nargs='*')
        add_dist_arg('-r', '--radius', help='Distribution description for the radius.')
        add_float_arg('-rho', '--density', help='Density of the particles.')
        add_float_arg('-mu', '--thermal-conductivity', help='Thermal conductivity of the particles.')

        add_float_arg('-bp', '--beam-power', help='Power of the laser beam in watts')
        add_float_arg('-br', '--beam-waist-radius', help='The waist radius of the LG01 beam (w_0)')

        add_float_arg('-gt', '--gas-temperature', help='The temperature of the surrounding gas [K].')
        add_float_arg('-gv', '--gas-viscosity', help='The dynamic viscosity of the surrounding gas [Pa*s].')
        add_float_arg('-gmu', '--gas-thermal-conductivity',
                      help='The thermal conductivity of the surrounding gas [W m^-1 K^-1].')
        add_float_arg('-gd', '--gas-density', help='The density of the surrounding gas  [kg m^-3].')

        add_arg('-G', '--gravity', help='Add a constant Z acceleration of -9.81 m/s^2.', action='store_true')
        add_float_arg('-b', '--boundary', help='Boundary (rmin, rmax, zmin, zmax) of the experiment.', nargs=4)

        parser.set_defaults(
            # Defaults given to me by Andrei Rode
            position=[{'kind': 'radial_gaussian', 'mu': 0.0, 'sigma': 20e-6}, 0.9e-3],
            velocity=[{'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.01}, {'kind': 'gaussian', 'mu': -10.0, 'sigma': 0.01}],
            radius={'kind': 'gaussian', 'mu': 2e-6, 'sigma': 10e-9},
            # Assume polystyrene spheres
            rho=1050.0,  # [kg/(m^3)]
            thermal_conductivity=0.030,  # [W/(m*K)]
            # Assume air at 293.15K, all values taken from https://www.engineeringtoolbox.com/. All values given are in
            # SI base units, e.g. the gas density is in [kg m^-3], the thermal conductivity in [W m^-1 K^-1)]
            gas_temperature=293.15, gas_viscosity=1.82e-5, gas_thermal_conductivity=25.87e-3, gas_density=1.204,
            gas_mass=4.81e-26,  # for an 'average' air particle (molar mass divided by Avogadro's constant)
            beam_power=1.0,  # [W]
            beam_waist_radius=4.0e-6,  # [m]
            boundary=[-1e-4, 1e-4, -1e-3, 1e-3]  # [m]
        )
        return parser

    @staticmethod
    def validate_args(args: argparse.Namespace):
        # Verify dimensionality match for position description and dimensions parameter
        if len(args.position) != 2:
            raise ValueError("The length of the position description must be 2!")
        if len(args.velocity) != 2:
            raise ValueError("The length of the velocity description must be 2!")

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
