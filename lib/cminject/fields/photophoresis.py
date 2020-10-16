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
Fields modeling photophoresis phenomena.
"""

from abc import ABC
from typing import Union, Tuple

import numba
import numpy as np
from scipy.constants import pi, Avogadro, R
from scipy.integrate import dblquad

from cminject.base import Field
from cminject.utils import empty_interval
from cminject.particles.spherical import ThermallyConductiveSphericalParticle
from .fluid_flow import DragForceInterpolationField


class VortexBeamPhotophoreticForceField(Field, ABC):
    """
    A base class that can calculate the irradiance of a Laguerre-Gaussian order 1 vortex beam, for a point in radial
    coordinates (r, z).
    """
    z_boundary = empty_interval

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
        z_0 = 2 * pi * self.beam_waist_radius**2 / self.beam_lambda
        # Note in Niko's thesis: "Change in waist over spherical surface may not be necessary"
        w = self.beam_waist_radius * np.sqrt(1 + (z / z_0)**2)
        intensity = (self.beam_power / pi) * (r**2 / w**4) * np.exp(-r**2 / w**2)
        return intensity


class DesyatnikovPhotophoreticLaserField(VortexBeamPhotophoreticForceField):
    """
    A Field that represents the acceleration exerted by an LG01 beam (in radial coordinates), according to the model
    by Desyatnikov, 2009.
    """
    def __init__(self,
                 gas_viscosity: float, gas_temperature: float, gas_thermal_conductivity: float,
                 gas_density: Union[float, Tuple[float, DragForceInterpolationField]], gas_mass: float,
                 beam_power: float, beam_waist_radius: float, beam_lambda: float = 523e-9,
                 z_position: float = 0.0):
        super().__init__(beam_power, beam_waist_radius, beam_lambda)
        self.gas_thermal_conductivity = gas_thermal_conductivity
        self.gas_density = gas_density
        self.gas_viscosity = gas_viscosity
        self.gas_temperature = gas_temperature
        self.molar_mass = gas_mass*Avogadro

        self.j1 = -0.5  # TODO assumed
        self.z0 = 2 * np.pi * self.beam_waist_radius ** 2 / self.beam_lambda  # see Desyatnikov 2009
        self.z_position = z_position or 0.0

    def kappa(self, particle_radius: float, particle_thermal_conductivity: float, gas_density: float):
        """The phenomenological constant \kappa, as described in Desyatnikov's 2009 paper."""
        return -self.j1 * 9*self.gas_viscosity**2 / \
               (2*particle_radius * gas_density * self.gas_temperature *
                (particle_thermal_conductivity + 2*self.gas_thermal_conductivity))

    def calculate_acceleration(self,
                               particle: ThermallyConductiveSphericalParticle,
                               time: float) -> np.array:
        a, r, z, m = particle.radius, particle.position[0], particle.position[1], particle.mass
        z = z - self.z_position  # For when the laser is offset in Z by some amount

        mu_a = particle.thermal_conductivity

        if type(self.gas_density) == float:
            gas_density = self.gas_density
        else:
            default_pressure, field = self.gas_density
            pressure = 0.0
            try:
                pressure = field.interpolate(particle.position)[2]  # pV=nRT
            except IndexError:
                pass
            pressure = pressure or default_pressure
            gas_density = self.molar_mass * pressure / (R * self.gas_temperature)

        f = self.pp_force(a, r, z, mu_a, gas_density)
        a = f / m

        return a

    @staticmethod
    @numba.jit("float64(float64,float64,float64,float64,float64)", nopython=True)
    def _transverse_integrand(y: float, x: float, a: float, r: float, w: float):
        return y / a * ((x**2 + (y + r)**2) / w**4) *\
               np.exp(-(x**2 + (y + r)**2) / w**2) *\
               a / np.sqrt(a**2 - x**2 - y**2)

    @staticmethod
    def _transverse(a: float, r: float, w: float):
        result, err = dblquad(
            DesyatnikovPhotophoreticLaserField._transverse_integrand,
            0, a,  # 0 to a since the integrand is even in x
            lambda x: -np.sqrt(a**2 - x**2), lambda x: np.sqrt(a**2 - x**2),
            args=(a, r, w),
        )
        return -1/np.pi * 2 * result  # 2*result since the integrand is even in x

    @staticmethod
    def _axial(a: float, w: float):
        faw = 1 - (1 + (a/w)**2) * np.exp(-(a/w)**2)  # = f(a**2/w**2), f taken from Desyatnikov's paper
        return faw

    def w(self, z: float):
        return self.beam_waist_radius * np.sqrt(1 + (z / self.z0) ** 2)

    def pp_force(self, a: float, r: float, z: float, k_f: float, rho: float) -> np.array:
        """
        Calculates the axial component of the photophoretic force, based on an approximation by Desyatnikov, 2009,
        for small particles with a << w.

        :param a: The particle radius.
        :param r: The particle's radial offset relative to the beam axis.
        :param z: The particle's axial offset relative to the beam origin.
        :param k_f: The particle's thermal conductivity.
        :param rho: The local density of the gas at (r, z).
        :return: A tuple containing the transverse (first element) and axial (second element) components of
            the photophoretic force.
        """
        P = self.beam_power
        kappa = self.kappa(a, k_f, rho)

        w = self.w(z)

        tr = self._transverse(a, r, w)
        ax = self._axial(a, w)

        return kappa * P * np.array([tr, ax])

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
