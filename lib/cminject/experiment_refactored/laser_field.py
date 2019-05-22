from functools import partial
from typing import Tuple

import numpy as np
from cminject.experiment_refactored.basic import SphericalParticle
from cminject.experiment_refactored.fluid_flow_field import FluidFlowField
from scipy.integrate import dblquad
from scipy.constants import R, pi

import cmath  # TODO can't we do this with np?

from cminject.experiment_refactored.base_classes import Field, Particle


class PhotophoreticLaserField(Field):
    def __init__(self, fluid_flow_field: FluidFlowField, power: float, lambda_: float, radius: float):
        self.fluid = fluid_flow_field
        self.power = power
        self.lambda_ = lambda_
        self.radius = radius

    def get_z_boundary(self) -> Tuple[float, float]:
        pass  # TODO

    def calculate_acceleration(self,
                               particle: SphericalParticle,
                               time: float) -> np.array:
        transverse, axial = self.calculate_photophoretic_force(particle)
        y_over_x = particle.position[1] / particle.position[0]
        return np.array([
            transverse * np.cos(np.arctan(y_over_x)),
            transverse * np.sin(np.arctan(y_over_x)),
            axial
        ])

    def calculate_vortex_intensity(self, x, y, z):
        """Vortex intensity distribution propagating along the z axis"""
        z_0 = pi * self.radius ** 2 / self.lambda_
        w = self.radius * np.sqrt(1 + (z / z_0))
        r2 = x ** 2 + y ** 2
        intensity = (4 * self.power / pi) * r2 / w ** 4 * np.exp(-2 * r2 / w ** 2)
        return intensity

    @staticmethod
    def calculate_sp_split(x, y, z, phi):
        """Returns the percentage of the S and P component in the beam."""
        r2 = x**2 + y**2 + z**2
        s = ((x * np.sin(phi) - y * np.cos(phi))**2) / r2
        p = ((x * np.cos(phi) - y * np.sin(phi))**2) / r2
        return s, p

    @staticmethod
    def calculate_absorption_coefficient(x, y, nt):
        """
        Calculates the absorption coefficients for S and P polarisation components according to the Fresnel equations.
        nt is the index of refraction of the sphere.
        """
        incident_angle = np.arcsin(np.sqrt(x**2 + y**2))
        if incident_angle != 0:
            transmit_angle = cmath.asin(np.sin(incident_angle) / nt)
            delta_angle = incident_angle - transmit_angle
            total_angle = incident_angle + transmit_angle
            p = 1 - abs(cmath.tan(delta_angle) / cmath.tan(total_angle))**2
            s = 1 - abs(cmath.sin(delta_angle) / cmath.sin(total_angle))**2
            return s, p, np.cos(incident_angle)
        else:
            p = 1 - abs((nt - 1) / (nt + 1))**2
            return p, p, 1.0

    def calculate_power_absorption(self, particle: Particle, x, y, z, nt):
        """
        The absorption is calculated for each component s and p.
        """
        p_abs, s_abs, cos_theta = self.calculate_absorption_coefficient(x, y, nt)

        p, s = self.calculate_sp_split(x, y, z, phi=0)
        p_abs *= p
        s_abs *= s

        offset_position = particle.position + np.array([x, y, z])
        intensity = self.calculate_vortex_intensity(*offset_position)
        absorbed_power = intensity * (p_abs + s_abs)
        return absorbed_power

    def calculate_integrated_absorption(self, particle: Particle, theta, phi, comp):
        """Integrate the absorption over a sphere"""
        x = np.sin(theta) * np.cos(phi) * self.radius
        y = np.sin(theta) * np.sin(phi) * self.radius
        z = np.cos(theta) * self.radius

        power = self.calculate_power_absorption(particle, x, y, z, 0)  # TODO nt?
        if comp == 'x':
            return abs(power * x) * np.sin(phi)
        elif comp == 'z':
            return power * z * np.sin(phi)  # TODO no abs?
        else:
            raise ValueError("comp needs to be 'x' or 'z'!")

    def calculate_photophoretic_force(self, particle: SphericalParticle):
        fl = self.fluid
        d = pi/2 * np.sqrt(pi/3) * fl.thermal_creep * (fl.dynamic_viscosity / fl.density) \
            * np.sqrt(fl.temperature * 8 * R / (pi * fl.molar_mass)) / fl.temperature
        char_p = (3/pi) * d * fl.temperature * particle.radius

        i = self.calculate_integrated_absorption()
        integrate_fn = partial(self.calculate_integrated_absorption, particle)
        first_half_sphere = dblquad(integrate_fn, 0, pi, 0, pi)[0]
        first_half_sphere_ax = dblquad(integrate_fn, 0, pi/2, 0, 2*pi)[0]
        second_half_sphere = dblquad(integrate_fn, 0, pi, pi, 2*pi)[0]
        second_half_sphere_ax = dblquad(integrate_fn, pi/2, pi, 0, 2 * pi)[0]

        diff_transverse = second_half_sphere - first_half_sphere
        diff_axial = second_half_sphere_ax - first_half_sphere_ax

        factor = d * particle.radius**2 * fl.pressure / (2 * char_p)  # TODO particle.kp ????
        return diff_transverse * factor, diff_axial * factor
