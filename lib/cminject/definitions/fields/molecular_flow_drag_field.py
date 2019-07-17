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

from typing import Tuple

import numpy as np

from scipy.constants import pi, Boltzmann
from scipy.special import erf,erfc

from cminject.definitions.base import Particle
from cminject.definitions.fields.regular_grid_interpolation_field import RegularGridInterpolationField
from cminject.definitions.particles import SphericalParticle


class MolecularFlowDragField(RegularGridInterpolationField):
    """
    A flow field that calculates a drag force exerted on a particle, based on the Epstein force for high velocities and interpolation
    on a grid defined by an HDF5 file like comsol_hdf5_tools.txt_to_hdf5 outputs.
    """
    def __init__(self, filename: str, m_gas: float = None, temperature: float = 4.0):
                 # For temperature calculations done by ParticleTemperaturePropertyUpdater
        # Store all the fixed initial properties
        self.temperature = temperature
        self.m_gas = m_gas

        # Short-term memoization storage
        self.interpolation_results = {}
        super().__init__(filename)


    def calculate_acceleration(self, particle: SphericalParticle, time: float) -> np.array:
        relative_velocity, _, pressure = self.interpolate(particle, time)
        f = self.calculate_drag_force(relative_velocity, pressure, particle)
        return f / particle.mass

    def is_particle_inside(self, particle: Particle, time: float) -> bool:
        return super().is_particle_inside(particle, time) and self.interpolate(particle, time)[2] > 0.0

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def interpolate(self, particle: Particle, time: float) -> Tuple[np.array, float, float]:
        """
        Calculates an (interpolated) relative velocity vector and pressure for a particle.
        Does memoization based on the particle identifier and current time.

        :param particle: The particle.
        :param time: The current time.
        :return: A 3-tuple of:
        - The velocity vector of the field relative to the particle,
        - The norm of this velocity vector,
        - The pressure exerted on the particle.
        """
        # Some basic memoizing so we can call this method often but have it calculated only once per time step
        # and particle.
        if particle.identifier in self.interpolation_results:
            stored_time, result = self.interpolation_results.get(particle.identifier)
            if time == stored_time:
                return result

        try:
            data = self.interpolator(tuple(particle.spatial_position,))
            pressure = data[self.number_of_dimensions]
            relative_velocity = data[:self.number_of_dimensions] - particle.velocity
            result = (relative_velocity, np.linalg.norm(relative_velocity), pressure)
        except ValueError:
            # Return zero velocity and zero pressure, so the particle will be considered outside
            result = np.zeros(self.number_of_dimensions), 0.0, 0.0

        # Store result for memoization and return it
        if particle.time_of_flight == time:
            self.interpolation_results[particle.identifier] = (particle.time_of_flight, result)
        return result

    def calculate_drag_force(self,
                             relative_velocity: np.array,
                             pressure: float,
                             particle: SphericalParticle) -> np.array:
        """
        Calculates the drag force using Epstein's law for spherical particles in molecular flow with corrections for high velocities

        :param relative_velocity: The velocity of the field relative to the particle
        :param pressure: The pressure exerted on the particle at the point
        :param particle: A Particle
        :return: The drag force as a (3,) np.array.
        """
        # Negative pressure or zero pressure is akin to the particle not being in the flow field, force is zero
        # and we can stop considering the particle to be inside this field
        if pressure <= 0:
            return np.zeros(self.number_of_dimensions)
        h=self.m_gas/(2*Boltzmann*self.temperature)
        # calculationg the force for two different cases 
        # 1. Epseints formula (V<10 m/s)
        # 2. Epsteins formula corrected for high velocities (V>=10 m/s)
        f_spec=np.where(abs(relative_velocity)<10,16/3*pressure*np.sqrt(pi*h)*particle.radius**2*relative_velocity,-pressure*np.sqrt(pi)*particle.radius**2/(2*h*relative_velocity**2)*(-2*np.exp(-h*relative_velocity**2)*np.sqrt(h)*relative_velocity*(1+2*h*relative_velocity**2)+np.sqrt(pi)*(1-4*h*relative_velocity**2-4*h**2*relative_velocity**4)*erf(np.sqrt(h)*relative_velocity)))
        force_vector = f_spec+1.8/3*pressure*pi**(3/2)*np.sqrt(h)*particle.radius**2*relative_velocity
        return force_vector


    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"This field was constructed from a {self.number_of_dimensions}-dimensional grid and is "
                f"thus incompatible with a {number_of_dimensions}-dimensional simulation."
            )


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
