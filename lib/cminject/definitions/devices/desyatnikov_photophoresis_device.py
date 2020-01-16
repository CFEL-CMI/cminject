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

from typing import List, Tuple, Union

import numpy as np
from cminject.definitions.property_updaters.brownian_motion import BrownianMotionPropertyUpdater
from scipy.constants import Boltzmann

from .base import Device
from cminject.definitions.fields.base import Field
from cminject.definitions.property_updaters.base import PropertyUpdater

from cminject.definitions.boundaries.cuboid import CuboidBoundary
from cminject.definitions.fields.fluid_flow import StokesDragForceField
from cminject.definitions.fields.photophoresis import DesyatnikovPhotophoreticLaserField
from cminject.definitions.fields.function_based import FunctionField

from cminject.definitions.particles.t_conductive_spherical import ThermallyConductiveSphericalParticle


# TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TODO this and the property updater should not live here! They should either not exist, or fall within
# TODO the family of drag force fields (just without grid interpolation). I chose to put them here for now, to allow
# TODO my combined simulations to finish faster without having to implement something new. This isn't good design,
# TODO it's only for me to be able to finish the thesis.     - Simon
# TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def get_uniform_drag_force(viscosity, velocity, temperature, m_gas, pressure):
    def uniform_drag_force(particle, time):
        knudsen = viscosity / (pressure * particle.radius) * \
                  np.sqrt(np.pi * Boltzmann * temperature / (2 * m_gas))
        s = 1 + knudsen * (1.231 + 0.4695 * np.exp(-1.1783 / knudsen))
        relative_velocity = velocity - particle.velocity

        force = 6 * np.pi * viscosity * particle.radius * relative_velocity * s
        return force / particle.mass
    return uniform_drag_force


class UniformBrownianMotionPropertyUpdater(PropertyUpdater):
    def __init__(self, viscosity, temperature, m_gas, pressure, dt):
        self.viscosity = viscosity
        self.temperature = temperature
        self.m_gas = m_gas
        self.pressure = pressure
        self.dt = dt

    def update(self, particle: ThermallyConductiveSphericalParticle, time) -> bool:
        if self.pressure >= 0.0:
            knudsen = self.viscosity / (self.pressure * particle.radius) * \
                      np.sqrt(np.pi * Boltzmann * self.temperature / (2 * self.m_gas))
            s = 1 + knudsen * (1.231 + 0.4695 * np.exp(-1.1783 / knudsen))

            s0 = 216 * self.viscosity * Boltzmann * self.temperature /\
                 (np.pi**2 * (2 * particle.radius)**5 * particle.rho**2 * s)
            a = BrownianMotionPropertyUpdater._generate_random_3d_vec() * np.sqrt(np.pi * s0 / self.dt)
            particle.position = BrownianMotionPropertyUpdater._generate_new_pos_2d(particle.position, a, self.dt)
            return True
        else:
            return False

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != 2:
            raise ValueError("Only available in n=2 dimensions!")


class DesyatnikovPhotophoresisDevice(Device):
    """
    A photophoretic LG01 vortex laser device based on Desyatnikov 2009.
    """
    def __init__(self, r_boundary: Tuple[float, float], z_boundary: Tuple[float, float],
                 gas_temperature: float, gas_viscosity: float, gas_thermal_conductivity: float,
                 gas_density: Union[float, Tuple[float, StokesDragForceField]], gas_mass: float,
                 beam_power: float, beam_waist_radius: float,
                 flow_gas_velocity: float = None, flow_gas_pressure: float = None, z_position: float = 0.0):
        """
        The constructor for DesyatnikovPhotophoresisDevice.

        :param r_boundary: The r interval where the laser is active
        :param z_boundary: The z interval where the laser is active
        :param gas_temperature: The temperature of the gas the laser is positioned in [K]
        :param gas_viscosity: The viscosity of the gas the laser is positioned in [Pa*s]
        :param gas_thermal_conductivity: The thermal conductivity of the gas the laser is positioned in [W/(m*K)]
        :param gas_density: The density of the gas the laser is positioned in, either as a fixed value or a tuple of a
            default value (a float) and a flow field to base the density on an interpolated pressure.
        :param gas_mass: The mass of a single gas molecule [kg]
        :param beam_power: The total power of the laser [W]
        :param beam_waist_radius: The radius of the laser beam waist [m]
        :param flow_gas_pressure: The pressure of a uniform drag force field
            (optional, but required if flow_gas_velocity is passed) [Pa]
        :param flow_gas_velocity: The velocity of a uniform drag force field
            (optional, but required if flow_gas_pressure is passed) [m/s]
        :param z_position: The Z position of the laser beam waist (0.0 by default)
        """

        if beam_power:
            pp_field = DesyatnikovPhotophoreticLaserField(
                gas_viscosity=gas_viscosity, gas_temperature=gas_temperature,
                gas_thermal_conductivity=gas_thermal_conductivity, gas_density=gas_density, gas_mass=gas_mass,
                beam_power=beam_power, beam_waist_radius=beam_waist_radius, z_position=z_position
            )
            fields: List[Field] = [pp_field]
        else:
            fields = []

        if flow_gas_velocity:
            drag_field = FunctionField(get_uniform_drag_force(gas_viscosity, np.array([0, flow_gas_velocity]),
                                                              gas_temperature, gas_mass, flow_gas_pressure))
            fields.append(drag_field)

        boundary: CuboidBoundary = CuboidBoundary(intervals=[r_boundary, z_boundary])
        super().__init__(fields=fields, boundary=boundary)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
