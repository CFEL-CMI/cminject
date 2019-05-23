"""
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
"""
from typing import Tuple

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.constants import pi, Avogadro, Boltzmann, R

from cminject.experiment_refactored.base_classes import Field, Particle
from cminject.experiment_refactored.basic import SphericalParticle, ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.comsol_hdf5_tools import hdf5_to_data_grid


class FluidFlowField(Field):
    """
    A flow field that calculates a drag force exerted on a particle, based on a grid defined by an HDF5 file
    that has been converted from a COMSOL .txt file, see comsol_hdf5_tools.txt_to_hdf5.
    """
    def __init__(self, filename: str,
                 density: float, dynamic_viscosity: float,
                 temperature: float = 4.0, pressure: float = 100.0, thermal_creep: float = 1.0,
                 # inflow_speed: float = 20.0, outflow_pressure: float = 0.0, k_n: float = 912.0,
                 # kinetic_d: float = 260.0e-12, conv: float = 0.00001,
                 # m_gas_mol: float = 0.004002602,
                 molar_mass: float = 0.004002602, molar_heat_capacity: float = 12.5,
                 specific_gas_constant: float = 2077.0, specific_heat_capacity: float = 3116.0,
                 m_gas: float = 6.6e-27, speed_of_sound: float = 117.7,
                 scale_slip: float = 1.0):
        # Store all the fixed initial properties
        self.density = density
        self.dynamic_viscosity = dynamic_viscosity
        self.temperature = temperature
        self.scale_slip = scale_slip
        self.pressure = pressure
        self.thermal_creep = thermal_creep
        self.molar_mass = molar_mass
        self.molar_heat_capacity = molar_heat_capacity
        self.specific_gas_constant = specific_gas_constant
        self.specific_heat_capacity = specific_heat_capacity
        self.m_gas = m_gas
        self.speed_of_sound = speed_of_sound

        # Storage to quickly remember and look up particles that are considered to be outside this field
        self.outside_particles = set()
        self.interpolation_results = {}

        # Construct the interpolator for the drag force from the HDF5 file passed by file name
        self.filename = filename
        x, y, z, data_grid = self._read_from_hdf5(self.filename)
        self.min_z, self.max_z = np.min(z), np.max(z)
        self.f_drag = RegularGridInterpolator((x, y, z), data_grid)

    @staticmethod
    def _read_from_hdf5(filename: str) -> Tuple[np.array, np.array, np.array, np.array]:
        print(f"Reading in HDF5 file {filename}...")
        x, y, z, data_grid = hdf5_to_data_grid(filename)
        print(f"Done reading in HDF5 file {filename}.")
        # Return x/y/z and the data grid. Both are needed to construct a RegularGridInterpolator.
        return x, y, z, data_grid

    def calculate_acceleration(self, particle: SphericalParticle, time: float) -> np.array:
        relative_velocity, _, pressure = self._interpolate(particle)
        return self.calculate_drag_force(relative_velocity, pressure, particle) / particle.mass

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.min_z, self.max_z

    def is_particle_inside(self, particle: Particle) -> bool:
        # NOTE/TODO: Particles are not able to leave the flow field and enter it again later!
        return self.min_z <= particle.position[2] <= self.max_z \
               and particle.identifier not in self.outside_particles

    def _interpolate(self, particle: Particle) -> Tuple[np.array, float, float]:
        """
        Calculates an (interpolated) relative velocity vector and pressure for a particle.
        Does memoization based on the particle identifier and current time.
        :param particle: The particle
        :return: The velocity vector of the field relative to the particle,
                 the norm of this vector, and the pressure exerted on the particle
        """
        # Some basic memoizing so we can call this method often but have it calculated only once per time step
        # and particle.
        if False and particle.identifier in self.interpolation_results:
            stored_time, result = self.interpolation_results.get(particle.identifier)
            if particle.time == stored_time:
                return result
            else:
                pass

        try:
            data = self.f_drag(tuple(particle.position))
            pressure = data[3]
            relative_velocity = data[:3] - particle.velocity
            result = (relative_velocity, np.linalg.norm(relative_velocity), pressure)
        except ValueError:
            # Return zero velocity and zero pressure, so the particle will be considered outside
            result = np.zeros(3), 0.0, 0.0

        # Store result for memoization and return it
        self.interpolation_results[particle.identifier] = (particle.time, result)
        return result

    def calculate_drag_force(self,
                             relative_velocity: np.array,
                             pressure: float,
                             particle: SphericalParticle) -> np.array:
        """
        Calculates the drag force using Stokes' law for spherical particles in continuum
        :param relative_velocity: The velocity of the field relative to the particle
        :param pressure:
        :param particle:
        :return:
        """
        # Negative pressure or zero pressure is akin to the particle not being in the flow field, force is zero
        # and we can stop considering the particle to be inside this field
        if pressure <= 0:
            self.outside_particles.add(particle.identifier)
            return np.zeros(3)

        force_vector = 6 * np.pi * self.dynamic_viscosity * particle.radius * relative_velocity
        return force_vector / self.calc_slip_correction(pressure=pressure, particle=particle)

    def calc_slip_correction(self, pressure: float, particle: SphericalParticle) -> float:
        """
        Calculates the slip correction factor with temperature corrections.
        The Sutherland constant for helium is 79.4 at reference temperature of 273.0 K.
        I took a reference pressure of 1 Pascal, in which the mean free path of helium is 0.01754.

        see J. Aerosol Sci. 1976 Vol. 7. pp 381-387 by Klaus Willeke
        """
        c_sutherland = 79.4
        ref_temperature = 273.0  # K

        factor = self.temperature / ref_temperature \
                 * (1 + c_sutherland / ref_temperature) \
                 / (1 + c_sutherland / self.temperature)
        knudsen = 0.01754 * 1 / (pressure * particle.radius) * factor

        s = 1 + knudsen * (1.246 + (0.42 * np.exp(-0.87 / knudsen)))
        return s * self.scale_slip

    def calc_mean_free_path(self):
        k = 1.38065e-23
        return self.dynamic_viscosity / self.pressure * np.sqrt(pi * k * self.temperature / (2 * self.m_gas))

    def calc_thermal_conductivity(self, particle: Particle) -> float:
        """
        Calculates the thermal conductivity based on the ideal gas law
        http://hyperphysics.phy-astr.gsu.edu/hbase/thermo/thercond.html
        :param particle: A Particle
        :return: The thermal conductivity
        """
        v: float = self._interpolate(particle)[1]
        mean_free_path = self.calc_mean_free_path()
        k = self.pressure * v * mean_free_path * self.molar_heat_capacity / \
            (3 * self.molar_mass * self.specific_gas_constant * self.temperature)
        return k

    def calc_prandtl(self, particle: Particle):
        """
        Calculates the Prandtl number.
        :param particle: A Particle
        :return: The Prandtl number.
        """
        return self.specific_heat_capacity * self.dynamic_viscosity / self.calc_thermal_conductivity(particle)

    def calc_reynolds(self, particle: Particle, d: float) -> float:
        """
        Calculates the Reynolds number.
        :param d: A characteristic length (in meters)
        :param particle: A Particle
        :return: The Reynolds number.
        """
        rho = self.pressure * self.specific_gas_constant * self.temperature
        v: float = self._interpolate(particle)[1]
        return rho * v * d / self.dynamic_viscosity

    def calc_nusselt(self, reynolds: float, prandtl: float, d: float, mu_s: float) -> float:
        """
        Calculates the Nusselt number based on the Reynolds and Prandtl number.
        :param d: A characteristic length (in meters)
        :param mu_s: The dynamic viscosity at the surface temperature of nanoparticle
        :param reynolds: The Reynolds number.
        :param prandtl: The Prandtl number.
        :return: The Nusselt number
        """
        re, pr = reynolds, prandtl  # to keep the equation shorter
        nu = 2 + (0.4 * re**0.5 + 0.06 * re**0.67) * pr**0.4 * (self.dynamic_viscosity / mu_s)**0.25
        return nu

    def calc_mach(self, particle: Particle):
        """
        Calculates the Mach number of the fluid at a Particle's position
        :param particle: A Particle
        :return: The Mach number
        """
        v: float = self._interpolate(particle)[1]
        return v / self.speed_of_sound

    def calc_convective_heat_transfer_coefficient(self, particle: Particle, nu: float, d: float):
        return nu / d * self.calc_thermal_conductivity(particle)

    def calc_particle_temperature(self, particle: ThermallyConductiveSphericalParticle, time: float):
        """
        Calculates a particle's "temperature" within this fluid at a given time.
        :param particle: A particle
        :param time: A time
        :return: The particle's temperature
        """
        d = 2 * particle.radius
        mu_s = 1.96e-5
        re = self.calc_reynolds(particle, d)
        pr = self.calc_prandtl(particle)
        m = self.calc_mach(particle)

        nu_0 = self.calc_nusselt(re, pr, d, mu_s)
        nu = nu_0 / (1 + 3.42 * (m * nu_0 / (re * pr)))

        h = self.calc_convective_heat_transfer_coefficient(particle, nu, d)
        c = particle.thermal_conductivity * particle.mass
        r = h * (4 * pi * particle.radius**2) / c  # TODO move surface area to particle class? method or field?
        return self.temperature + (particle.temperature - self.temperature) * np.exp(-r * time)

    def calc_particle_collisions(self,
                                 particle: ThermallyConductiveSphericalParticle,
                                 dt: float) -> Tuple[float, float]:
        """
        Estimates the number of collisions of a particle within a time frame based on Kinetic Gas Theory
        :param particle: A Particle
        :param dt: A time delta
        :return: A tuple of (the number of collisions, the new temperature based on them)
        """
        v = np.sqrt(8 * self.temperature * Boltzmann / (self.m_gas * pi))
        v += self._interpolate(particle)[1]
        n = self.pressure * Avogadro / (R * self.temperature)
        c = 0.5 * n * v * dt * (4 * pi * particle.radius**2)
        k = (particle.mass + self.m_gas)**2 / (2 * particle.mass * self.m_gas)
        # TODO in the original code, dT was calculated but unused?
        collision_temperature = self.temperature + (particle.temperature - self.temperature) * np.exp(-c / k)
        return c, collision_temperature

"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
