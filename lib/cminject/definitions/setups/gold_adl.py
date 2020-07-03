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
from typing import Tuple

import numpy as np

from cminject.definitions.boundaries.grid_field_based import GridFieldBasedBoundary
from cminject.definitions.detectors.simple_z import SimpleZDetector
from cminject.definitions.devices.base import Device, Boundary
from cminject.definitions.devices.desyatnikov_photophoresis_device import DesyatnikovPhotophoresisDevice
from cminject.definitions.fields.fluid_flow import StokesDragForceField
from cminject.definitions.particles.t_conductive_spherical import ThermallyConductiveSphericalParticle
from cminject.definitions.property_updaters.brownian_motion import BrownianMotionPropertyUpdater
from cminject.definitions.setups import Setup
from cminject.definitions.sources.variable_distributions import VariableDistributionSource
from cminject.experiment import Experiment
from cminject.utils.args import SetupArgumentParser


class SkimmersBoundary(Boundary):
    def __init__(self, first_skimmer_flow_field: StokesDragForceField,
                 inter_distance: float, skimmer_length: float, tube_length: float,
                 skimmer_min_radius: float, skimmer_max_radius: float,
                 skimmer_min_z: float):
        self.first_skimmer_flow_field = first_skimmer_flow_field
        self.inter_distance = inter_distance
        self.skimmer_length = skimmer_length
        self.tube_length = tube_length
        self.skimmer_min_radius = skimmer_min_radius
        self.skimmer_max_radius = skimmer_max_radius

        self._z_boundary = (skimmer_min_z, skimmer_length + tube_length)

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        r, z = position
        r = np.abs(r)

        zmin, zmax = self.z_boundary
        if not (zmin <= z <= zmax):
            return False

        if z <= 0.0:
            return self.first_skimmer_flow_field.is_particle_inside(position, time)
        elif z <= self.inter_distance:
            return True  # In the region in between we never consider particles to be lost
        elif self.inter_distance < z <= self.skimmer_length:
            # In the trapezoid region, verify the particle is within the offset triangle
            lng = self.skimmer_length
            r0, r1 = self.skimmer_min_radius, self.skimmer_max_radius
            return z >= lng*((r-r0) / (r1-r0))
        else:
            # In the rest of the Z boundary of this device, just require that the particles are inside the tube
            return r <= self.skimmer_max_radius

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary


NITROGEN_PARAMS = {
    'dynamic_viscosity': 1.76e-5,  # [Pa*s]
    'm_gas': 4.7e-26,  # [kg]
    'temperature': 293.15,  # [K]
    'slip_correction_model': 'room_temp'
}


class SkimmersDevice(Device):
    def __init__(self, filename, inter_distance: float, skimmer_length: float, tube_length: float,
                 skimmer_min_r: float, skimmer_max_r: float, skimmer_min_z: float):
        field = StokesDragForceField(filename=filename, **NITROGEN_PARAMS)
        boundary = SkimmersBoundary(
            first_skimmer_flow_field=field,
            inter_distance=inter_distance,
            skimmer_length=skimmer_length,
            tube_length=tube_length,
            skimmer_min_radius=skimmer_min_r,
            skimmer_max_radius=skimmer_max_r,
            skimmer_min_z=skimmer_min_z
        )

        super().__init__(fields=[field], boundary=boundary)


class ADLStackDevice(Device):
    def __init__(self, filename, z_offset):
        field = StokesDragForceField(filename=filename, **NITROGEN_PARAMS, offset=np.array([0, z_offset]))
        boundary = GridFieldBasedBoundary(field)
        super().__init__(fields=[field], boundary=boundary)


class GoldADLSetup(Setup):
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        skimmer_flow_file = args.skimmer1_file
        adl_flow_file = args.adl_file

        inter_distance = 2.95e-3
        skimmer_length = 2e-2  # TODO don't know if this is right ???
        skimmer_tube_length = 15e-2
        skimmer_min_z = -2.35835e-2
        skimmer_min_r = 0.5e-3
        skimmer_max_r = 2.0e-2

        exit_offset = 0.0003
        adl_extent = 0.4272
        wtf_offset = -1e-3  # TODO i know this is necessary, but where is it coming from?
        adl_offset = exit_offset + adl_extent + wtf_offset + skimmer_length + inter_distance + skimmer_tube_length

        adl_exit_position = 0.59915  # 0.60015 is where the entire ADL flow field terminates

        devices = [
            SkimmersDevice(
                skimmer_flow_file, inter_distance=inter_distance,
                skimmer_length=skimmer_length, tube_length=skimmer_tube_length,
                skimmer_min_r=skimmer_min_r, skimmer_max_r=skimmer_max_r, skimmer_min_z=skimmer_min_z
            ),
            ADLStackDevice(adl_flow_file, z_offset=adl_offset),
        ]

        if args.beam_power is not None:  # if missing or 0.0, we don't need to create a device
            devices += [DesyatnikovPhotophoresisDevice(
                # All gas properties are assuming dilute nitrogen at 293.15K.
                gas_density=(0.059, devices[1].fields[0]),
                # taken for dilute nitrogen from Lemmon and Jacobsen 2003
                gas_viscosity=17.8771e-6,  # [Pa*s]
                gas_thermal_conductivity=25.9361e-3,  # [W/(m*K)]

                gas_mass=4.652e-26,  # [kg], one molecule
                gas_temperature=293.15,  # [K]

                beam_waist_radius=args.beam_waist_radius,  # [m]
                beam_power=args.beam_power,  # [W]
                r_boundary=(-1e-3, 1e-3),  # [m]
                z_boundary=(adl_exit_position, adl_exit_position + 1e-1),  # [m]
                z_position=adl_exit_position + 1e-2,  # [m]
            )]

        detectors = [SimpleZDetector(i, adl_exit_position + i*5e-4) for i in range(20)]
        sources = [
            VariableDistributionSource(
                subclass=ThermallyConductiveSphericalParticle,
                number_of_particles=main_args.nof_particles,
                rho=19320.0,  # Assuming 50nm gold particles
                radius=50e-9,
                thermal_conductivity=315.0,  # [W / (m*K)], taken from Wikipedia
                specific_heat=0.0,
                temperature=293.15,
                #rho=1050.0,  # Assuming expanded Polystyrene at 50nm
                #radius=50e-9,
                #thermal_conductivity=0.030,

                position=[{'kind': 'radial_gaussian', 'mu': 0.0, 'sigma': 3.0e-3},
                          skimmer_min_z+abs(skimmer_min_z*0.001)],
                velocity=[{'kind': 'gaussian', 'mu': 1e-3, 'sigma': 1e-5},
                          {'kind': 'gaussian', 'mu': 0.155, 'sigma': 0.001}]
                # for ADL exit, approximated:
                # position=[{'kind': 'gaussian', 'mu': 0.0, 'sigma': 10e-6}, 0.6],
                # velocity=[{'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.1}, {'kind': 'gaussian', 'mu': 39.0, 'sigma': 0.01}]
            )
        ]
        # Construct the Experiment, still missing brownian motion
        exp = Experiment(
            detectors=detectors,
            sources=sources,
            devices=devices,
            time_step=1e-5,
            number_of_dimensions=2
        )

        if args.brownian:
            # Add brownian motion property updaters _after_ the Experiment has decided on a time step to use.
            # TODO Might be a good idea to refactor this automatic choice of time step to somewhere else.
            dt = exp.time_step
            bm_updaters = [BrownianMotionPropertyUpdater(field=devices[0].fields[0], dt=dt),
                           BrownianMotionPropertyUpdater(field=devices[1].fields[0], dt=dt)]
            exp.property_updaters = bm_updaters + exp.property_updaters

        return exp

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()

        parser.add_argument('-fS1', '--skimmer1-file', help='The flow field file (HDF5) for the first skimmer',
                            required=True, type=str)
        parser.add_argument('-fADL', '--adl-file', help='The flow field file (HDF5) for the ADL setup',
                            required=True, type=str)

        parser.add_argument('-bP', '--beam-power', help='Beam power of a photophoretic LG01 vortex laser beam [W]',
                            type=float)
        parser.add_argument('-bw0', '--beam-waist-radius',
                            help='Beam waist radius of a photophoretic LG01 vortex laser beam [m]',
                            type=float, default=8e-6)

        parser.add_argument('-B', '--brownian', help='Enable brownian motion', action='store_true')

        return parser

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
