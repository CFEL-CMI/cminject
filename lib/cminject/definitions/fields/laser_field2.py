from typing import Tuple

import numpy as np
from cminject.definitions import Field
from cminject.definitions.particles import ThermallyConductiveSphericalParticle


class PhotophoreticLaserField2(Field):
    def __init__(self, w0: float, power: float, radius: float, position: np.array, lambda_: float = 532e-9):
        self.w0 = w0
        self.power = power
        self.number_of_dimensions = 3
        self.position = position
        self.radius = radius
        self.lambda_ = lambda_

    def calculate_acceleration(self, particle: ThermallyConductiveSphericalParticle, time: float) -> np.array:
        return np.zeros(self.number_of_dimensions)

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return (0, float('inf'))

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != 3:
            raise ValueError(f"Can only simulate {self.__class__} in 3D!")
        # TODO allow in 2D too

    def vortex_intensity(self, particle: ThermallyConductiveSphericalParticle) -> Tuple[float, float]:
        """
        Calculates the vortex intensity (propagating along the Z axis) at a particle's x/y/z position
        :param particle: The particle whose position to look at
        :return: A 2-tuple of the intensity and the power on the particle
        """
        x, y, z = particle.spatial_position
        w0, lambda_, power = self.w0, self.lambda_, self.power
        z0 = np.pi * w0**2 / lambda_
        w = w0 * np.sqrt(1 + (z / z0)**2)

        particle_power = power * (1 - (1+(1/w0)**2) * np.exp(-1 / (w0**2)))

        r_squared = x**2 + y**2
        intensity = (4 * power / np.pi) * r_squared / (w**4) * np.exp(-2 * r_squared / (w**2))
        intensity = intensity / particle.radius**2
        return intensity, particle_power

    def component_polarizations(self, x: float, y: float, phi: float) -> Tuple[float, float]:
        """
        Calculates the polarization components on a sphere based on an x/y position on the sphere
        and the angle between phase components of the incoming beam
        :param x: The x position on the sphere
        :param y: The y position on the sphere
        :param phi: The angle between phase components of the incoming beam
        :return: A 2-tuple of the p- and s-component
        """
        sum = x**2 + y**2
        if sum > 1:
            p_comp, s_comp = 0, 0
        elif sum == 0:
            p_comp, s_comp = 1, 0
        else:
            p_comp = ((x * np.cos(phi) - y * np.sin(phi)) ** 2) / sum
            s_comp = ((x * np.sin(phi) - y * np.cos(phi))**2) / sum

        return p_comp, s_comp

    def component_absorptions(self, x: float, y: float, nt: float) -> Tuple[float, float]:
        """
        Calculates the absorption coefficients for S and P polarisation components according to the
        Fresnel formulae.
        :param x: The x position on the sphere
        :param y: The y position on the sphere
        :param nt: The refractive index of the sphere
        :return: A 2-tuple of the absorption of the p- and s-component
        """
        incident_angle = np.arcsin(np.sqrt(x**2 + y**2))
        transmit_angle = np.arcsin(np.sin(incident_angle / nt))

        diff = incident_angle - transmit_angle
        sum = incident_angle + transmit_angle

        if incident_angle == 0:
            p = 1 - np.abs((nt - 1)/(nt + 1))**2
            s = p
        else:
            p = 1 - np.abs(np.tan(diff) / np.tan(sum)) ** 2
            s = 1 - np.abs(np.sin(diff) / np.sin(sum)) ** 2

        return p, s

    def power_absorption(self, particle: ThermallyConductiveSphericalParticle, phi: float, alpha: float = 1.0):
        """
        Calculates the power absorbed by a particle at some position for a given
        angle between linear polarisation states.
        :param particle: The particle to calculate the absorbed power for.
        :param phi: The angle between linear polarisation states
        :param alpha: The absorption coefficient of the particle. 1.0 by default.
        :return: TODO
        """
        x, y, z = particle.spatial_position
        p, s = self.component_polarizations(x, y, phi)
        p_abs, s_abs = self.component_absorptions(x, y, particle.refractive_index)
        s_absorption_strength = s * s_abs * z
        p_absorption_strength = p * p_abs * z

        # Scaling to W / cm^2
        intensity = self.vortex_intensity(particle) / 1e4
        # TODO ... if z <