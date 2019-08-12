#!/usr/bin/env python
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

from abc import ABC

import numpy as np
from cminject.definitions.base import Field, empty_interval
from cminject.definitions.particles import ThermallyConductiveSphericalParticle
from scipy.constants import pi


class VortexBeamPhotophoreticForceField(Field, ABC):
    z_boundary = empty_interval

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != 2:
            raise ValueError("Photophoretic force fields can currently only handle a 2D simulation setup.")

    def __init__(self, beam_power: float, beam_waist_radius: float, beam_lambda: float = 523e-9):
        self.beam_power = beam_power
        self.beam_lambda = beam_lambda
        self.beam_waist_radius = beam_waist_radius

    def calculate_vortex_irradiance(self, r: float, z: float):
        """
        Vortex irradiance distribution for a vortex beam propagating from (0,0) along the z axis towards the positive.
        The beam is radially symmetric around the Z axis, so only r and z are passed.

        :param r: The polar radius in the plane transverse to the optical Z axis
        :param z: The position on the optical Z axis
        """
        z_0 = pi * self.beam_waist_radius ** 2 / self.beam_lambda
        # Note in Niko's thesis: "Change in waist over spherical surface may not be necessary"
        w = self.beam_waist_radius * np.sqrt(1 + (z / z_0) ** 2)
        intensity = (4 * self.beam_power / pi) * r**2 / w**4 * np.exp(-2 * r**2 / w**2)
        return intensity


class DesyatnikovPhotophoreticLaserField(VortexBeamPhotophoreticForceField):
    def __init__(self,
                 gas_viscosity: float, gas_temperature: float,
                 gas_thermal_conductivity: float, gas_density: float,
                 beam_power: float, beam_waist_radius: float, beam_lambda: float = 523e-9):
        super().__init__(beam_power, beam_waist_radius, beam_lambda)
        self.gas_thermal_conductivity = gas_thermal_conductivity
        self.gas_density = gas_density
        self.gas_viscosity = gas_viscosity
        self.gas_temperature = gas_temperature

        self.j1 = -0.5  # TODO assumed
        self.z0 = 2 * np.pi * self.beam_waist_radius ** 2 / self.beam_lambda  # see Desyatnikov 2009

    def kappa(self, particle_radius: float, particle_thermal_conductivity: float):
        return -self.j1 * 9*self.gas_viscosity**2 / \
               (2*particle_radius * self.gas_density * self.gas_temperature *
                (particle_thermal_conductivity + 2*self.gas_thermal_conductivity))

    def calculate_acceleration(self,
                               particle: ThermallyConductiveSphericalParticle,
                               time: float) -> np.array:
        a, r, z, m = particle.radius, particle.position[0], particle.position[1], particle.mass
        mu_a = particle.thermal_conductivity
        return self.pp_force(a, r, z, mu_a) / m

    def pp_force(self, a: float, r: float, z: float, mu_a: float) -> np.array:
        """
        Calculates the axial component of the photophoretic force, based on an approximation by Desyatnikov, 2009,
        for small particles with a << w.

        :param a: The particle radius.
        :param r: The particle's radial offset relative to the beam axis.
        :param z: The particle's axial offset relative to the beam origin.
        :param mu_a: The particle's thermal conductivity.
        :return: A tuple containing the transverse (first element) and axial (second element) components of
            the photophoretic force.
        """
        P = self.beam_power
        kappa = self.kappa(a, mu_a)
        w = self.beam_waist_radius * np.sqrt(1 + (z / self.z0)**4)
        return np.array([
            (-kappa * P * 4/3 * r * a**3) / w**4,
            kappa * P/2 * (a / w)**4
        ])


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""