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

from cminject.base import Particle, Action, Setup
from cminject.detectors.simple_z import SimpleZDetector
from cminject.devices.desyatnikov_photophoresis_device import DesyatnikovPhotophoresisDevice
from cminject.particles.spherical import ThermallyConductiveSphericalParticle
from cminject.sources.variable_distributions import VariableDistributionSource
from cminject.experiment import Experiment
from cminject.utils.args import dist_description, SetupArgumentParser

M_GAS_AIR = 5.6e-26


class MarkAsLostWhenVZIsPositive(Action):
    """
    Marks particles as lost when their z coordinate and v_z are both positive.
    Since the beam only drives particles away from the source position in this area, we can consider the particles lost.
    """
    def __call__(self, particle: Particle, time: float) -> bool:
        if particle.position[-1] > 0.0 and particle.velocity[-1] > 0.0:
            particle.lost = True
        return False


class DesyatnikovPhotophoresisSetup(Setup):
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        experiment = Experiment(
            time_interval=main_args.time_interval, time_step=main_args.time_step,
            random_seed=main_args.seed, number_of_dimensions=2
        )

        experiment.add_source(VariableDistributionSource(
            main_args.nof_particles, position=args.position, velocity=args.velocity, radius=args.radius, density=args.rho,
            subclass=ThermallyConductiveSphericalParticle,
            subclass_kwargs={
                'specific_heat': 0.0, 'temperature': 293.15, 'thermal_conductivity': args.thermal_conductivity
            }
        ))

        experiment.add_device(DesyatnikovPhotophoresisDevice(
            r_boundary=args.boundary[:2], z_boundary=args.boundary[2:],
            gas_temperature=args.gas_temperature, gas_viscosity=args.gas_viscosity,
            gas_thermal_conductivity=args.gas_thermal_conductivity, gas_density=args.gas_density, gas_mass=M_GAS_AIR,
            beam_power=args.beam_power, beam_waist_radius=args.beam_waist_radius,
        ))

        for i, pos in enumerate(args.detectors):
            experiment.add_detector(SimpleZDetector(identifier=i, z_position=pos))

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
        parser.add_argument('-r', '--radius', help='Distribution description for the radius.', type=dist_description)
        parser.add_argument('-rho', '--density', help='Density of the particles.', type=float)
        parser.add_argument('-mu', '--thermal-conductivity', help='Thermal conductivity of the particles.', type=float)

        parser.add_argument('-bp', '--beam-power', help='Power of the laser beam in watts', type=float)
        parser.add_argument('-br', '--beam-waist-radius', help='The waist radius of the LG01 beam (w_0)', type=float)

        parser.add_argument('-gt', '--gas-temperature', help='The temperature of the surrounding gas [K].', type=float)
        parser.add_argument('-gv', '--gas-viscosity', help='The dynamic viscosity of the surrounding gas [Pa*s].',
                            type=float)
        parser.add_argument('-gmu', '--gas-thermal-conductivity',
                            help='The thermal conductivity of the surrounding gas [W m^-1 K^-1].', type=float)
        parser.add_argument('-gd', '--gas-density', help='The density of the surrounding gas  [kg m^-3].', type=float)

        parser.add_argument('-b', '--boundary', help='Boundary (rmin, rmax, zmin, zmax) of the experiment.',
                            type=float, nargs=4)

        parser.set_defaults(
            # Defaults given to me by Andrei Rode
            position=[{'kind': 'radial_gaussian', 'mu': 0.0, 'sigma': 20e-6},#42.47e-6},
                      0.9e-3],
            velocity=[{'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.01},
                      {'kind': 'gaussian', 'mu': -10.0, 'sigma': 0.01}],
            radius={'kind': 'gaussian', 'mu': 2e-6, 'sigma': 10e-9},

            # Add a uniform flow field
            uniform_flow_velocity=None,

            # Assume polystyrene spheres
            rho=1050.0,  # [kg/(m^3)]
            thermal_conductivity=0.030,  # [W/(m*K)]

            # Assume air at 293.15K, all values taken from https://www.engineeringtoolbox.com/
            gas_temperature=293.15,  # [K]
            gas_viscosity=1.82e-5,  # [Pa*s]
            gas_thermal_conductivity=25.87e-3,  # [W/(m*K)]
            gas_density=1.204,  # [kg/(m^3)]

            beam_power=1.0,  # [W]
            beam_waist_radius=4.0e-6,  # [m]

            boundary=[-1e-4, 1e-4, -1e-3, 1e-3],  # [m]
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
