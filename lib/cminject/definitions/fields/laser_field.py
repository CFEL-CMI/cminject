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

from functools import partial
from typing import Tuple

import numpy as np
from scipy.integrate import dblquad
from scipy.constants import R, pi

import cmath  # TODO can't we do this with np?

from cminject.definitions.base import Field, empty_interval
from cminject.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.definitions.fields.stokes_fluid_flow_field import StokesFluidFlowField


class PhotophoreticLaserField(Field):
    def __init__(self, fluid_flow_field: StokesFluidFlowField, power: float, radius: float, lambda_: float = 523e-9):
        self.fluid = fluid_flow_field
        self.power = power
        self.lambda_ = lambda_
        self.radius = radius

    def calculate_acceleration(self,
                               particle: ThermallyConductiveSphericalParticle,
                               time: float) -> np.array:
        transverse, axial = self.calculate_photophoretic_force(particle)
        y_over_x = particle.position[1] / particle.position[0]
        return np.array([
            transverse * np.cos(np.arctan(y_over_x)),
            transverse * np.sin(np.arctan(y_over_x)),
            axial
        ]) / particle.mass

    def calculate_vortex_intensity(self, x, y, z):
        """Vortex intensity distribution propagating along the z axis"""
        z_0 = pi * self.radius ** 2 / self.lambda_
        w = self.radius * np.sqrt(1 + (z / z_0))
        r2 = x ** 2 + y ** 2
        intensity = (4 * self.power / pi) * r2 / w ** 4 * np.exp(-2 * r2 / w ** 2)
        return intensity

    @staticmethod
    def calculate_sp_split(x, y, z, phi):
        """Returns the percentages of the S and P component in the beam."""
        r2 = x**2 + y**2 + z**2  # FIXME really z too? it's not like this in Niko's MATLAB code
        s = ((x * np.sin(phi) - y * np.cos(phi))**2) / r2
        p = ((x * np.cos(phi) - y * np.sin(phi))**2) / r2
        return s, p

    @staticmethod
    def calculate_absorption_coefficient(x: float, y: float, refraction_index: float = 1.0) \
            -> Tuple[float, float, float]:
        """
        Calculates the absorption coefficients for S and P polarisation components according to the Fresnel equations.
        :param x: The x position.
        :param y: The y position.
        :param refraction_index: The refraction index of the sphere.
        :return: A 3-tuple of the s and p absorption coefficients and the cos of the incident angle
        """
        incident_angle = np.arcsin(np.sqrt(x**2 + y**2))
        if incident_angle != 0.0:
            transmit_angle = cmath.asin(np.sin(incident_angle) / refraction_index)
            delta_angle = incident_angle - transmit_angle
            total_angle = incident_angle + transmit_angle
            p = 1 - abs(cmath.tan(delta_angle) / cmath.tan(total_angle))**2
            s = 1 - abs(cmath.sin(delta_angle) / cmath.sin(total_angle))**2
            return s, p, np.cos(incident_angle)
        else:
            p = 1 - abs((refraction_index - 1) / (refraction_index + 1))**2
            s = p
            return s, p, 1.0

    def calculate_power_absorption(self, particle: ThermallyConductiveSphericalParticle, x, y, z):
        """
        The absorption is calculated for each component s and p.
        """
        p_abs, s_abs, cos_theta = self.calculate_absorption_coefficient(x, y, particle.refractive_index)

        p, s = self.calculate_sp_split(x, y, z, phi=0)
        p_abs *= p
        s_abs *= s

        offset_position = particle.spatial_position + np.array([x, y, z])
        intensity = self.calculate_vortex_intensity(*offset_position) / self.radius**2
        absorbed_power = intensity * (p_abs + s_abs)
        return absorbed_power

    def absorption_integrator_fn(self, particle: ThermallyConductiveSphericalParticle, theta, phi, comp):
        """Integrate the absorption over a sphere"""
        x = np.sin(theta) * np.cos(phi) * particle.radius
        y = np.sin(theta) * np.sin(phi) * particle.radius
        z = np.cos(theta) * particle.radius
        power = self.calculate_power_absorption(particle, x, y, z)

        if comp == 'x':
            return abs(power * x) * np.sin(phi)
        elif comp == 'z':
            return abs(power * z) * np.sin(phi)
        else:
            raise ValueError("comp needs to be 'x' or 'z'!")

    def intensity_integrator_fn(self, particle: ThermallyConductiveSphericalParticle, theta, phi):
        """
        The function to be integrated to calculate the total intensity across a sphere section.
        :param particle: The spherical particle to integrate over.
        :param theta: Theta in spherical coordinates.
        :param phi: Phi in spherical coordinates.
        :return: The local intensity.
        """
        x = np.sin(theta) * np.cos(phi) * particle.radius
        y = np.sin(theta) * np.sin(phi) * particle.radius
        z = np.cos(theta)
        return np.sin(phi) * self.calculate_vortex_intensity(
            *(particle.position + np.array([x, y, z]))
        )

    def calculate_photophoretic_force(self, particle: ThermallyConductiveSphericalParticle) -> Tuple[float, float]:
        """
        Calculates the photophoretic force exerted on a particle based on a simple model,
        based amongst other things on the fluid flow force field this field was constructed with.
        :param particle: The particle instance.
        :return: The photophoretic forces in transverse and axial direction.
        """
        fl = self.fluid
        d = pi/2 * np.sqrt(pi/3) * fl.thermal_creep * (fl.dynamic_viscosity / fl.density) \
            * np.sqrt(fl.temperature * 8 * R / (pi * fl.molar_mass)) / fl.temperature
        char_p = (3/pi) * d * fl.temperature * particle.radius

        integrator_fn_ax = partial(self.absorption_integrator_fn, particle, comp='z')
        integrator_fn_tr = partial(self.absorption_integrator_fn, particle, comp='x')
        first_half_sphere = dblquad(integrator_fn_tr,     0,     pi,    0,   pi)[0]
        first_half_sphere_ax = dblquad(integrator_fn_ax,  0,     pi/2,  0,   2*pi)[0]
        second_half_sphere = dblquad(integrator_fn_tr,    0,     pi,    pi,  2*pi)[0]
        second_half_sphere_ax = dblquad(integrator_fn_ax, pi/2,  pi,    0,   2*pi)[0]

        diff_transverse = second_half_sphere - first_half_sphere
        diff_axial = second_half_sphere_ax - first_half_sphere_ax

        factor = d * particle.radius**2 * fl.pressure / (2 * particle.thermal_conductivity * char_p)
        return diff_transverse * factor, diff_axial * factor

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return empty_interval  # TODO this ain't great

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != 3:
            raise ValueError("PhotophoreticLaserField can currently only handle a 3D simulation setup.")


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""