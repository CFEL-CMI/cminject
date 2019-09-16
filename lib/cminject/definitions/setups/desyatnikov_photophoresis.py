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

from cminject.definitions import PropertyUpdater, Particle
from cminject.definitions.setups import Setup
from cminject.definitions.detectors import SimpleZDetector
from cminject.definitions.devices.desyatnikov_photophoresis_device import DesyatnikovPhotophoresisDevice, \
    UniformBrownianMotionPropertyUpdater
from cminject.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.definitions.sources import VariableDistributionSource
from cminject.experiment import Experiment
from cminject.utils.args import dist_description, SetupArgumentParser, auto_time_step


M_GAS_AIR = 5.6e-26


class MarkAsLostWhenVZIsPositivePropertyUpdater(PropertyUpdater):
    """
    Marks particles as lost when their z coordinate and v_z are both positive.
    Since the beam only drives particles away from the source position in this area, we can consider the particles lost.
    """
    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass

    def update(self, particle: Particle, time: float) -> bool:
        if particle.spatial_position[-1] > 0.0 and particle.velocity[-1] > 0.0:
            particle.lost = True
        return False


class DesyatnikovPhotophoresisSetup(Setup):
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        dt = main_args.time_step or auto_time_step(abs(args.velocity[-1]['mu']))

        devices = [DesyatnikovPhotophoresisDevice(
            r_boundary=args.boundary[:2], z_boundary=args.boundary[2:],
            gas_temperature=args.gas_temperature, gas_viscosity=args.gas_viscosity, flow_gas_pressure=args.gas_pressure,
            gas_thermal_conductivity=args.gas_thermal_conductivity, gas_density=args.gas_density, gas_mass=M_GAS_AIR,
            beam_power=args.beam_power, beam_waist_radius=args.beam_waist_radius,
            flow_gas_velocity=args.uniform_flow_velocity
        )]
        detectors = [SimpleZDetector(identifier=i, z_position=pos) for i, pos in enumerate(args.detectors)]
        sources = [VariableDistributionSource(
            main_args.nof_particles, position=args.position, velocity=args.velocity, radius=args.radius, rho=args.rho,
            subclass=ThermallyConductiveSphericalParticle, specific_heat=0.0, temperature=293.15,
            thermal_conductivity=args.thermal_conductivity
        )]

        if args.brownian:
            property_updaters = [
                # MarkAsLostWhenVZIsPositivePropertyUpdater(),
                UniformBrownianMotionPropertyUpdater(
                    viscosity=args.gas_viscosity, temperature=args.gas_temperature,
                    m_gas=M_GAS_AIR, pressure=args.gas_pressure, dt=dt
                ),
            ]
        else:
            property_updaters = []

        return Experiment(
            devices=devices, detectors=detectors, sources=sources, property_updaters=property_updaters,
            time_interval=main_args.time_interval, time_step=dt,
            delta_z_end=0.0, seed=main_args.seed, number_of_dimensions=2
        )

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
        parser.add_argument('-gp', '--gas-pressure', help='The pressure of the surrounding gas [Pa].', type=float)

        parser.add_argument('-F', '--uniform-flow-velocity', help='Add a uniform flow field with a given Z velocity',
                            type=float, required=False)
        parser.add_argument('-B', '--brownian',
                            help='Enable brownian motion in a uniform drag force flow field. Requires -F to be set.',
                            action='store_true')

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
            gas_pressure=10000.0,  # [Pa]

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
