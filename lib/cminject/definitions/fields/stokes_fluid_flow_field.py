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

from cminject.definitions.base import Particle
from cminject.definitions.fields.regular_grid_interpolation_field import RegularGridInterpolationField
from cminject.definitions.particles import SphericalParticle


class StokesFluidFlowField(RegularGridInterpolationField):
    """
    A flow field that calculates a drag force exerted on a particle, based on the stokes force and interpolation
    on a grid defined by an HDF5 file like comsol_hdf5_tools.txt_to_hdf5 outputs.
    """
    def __init__(self, filename: str,
                 dynamic_viscosity: float, density: float, m_gas: float = None, temperature: float = 4.0,
                 slip_correction_model: str = None, slip_correction_scale: float = 1.0,
                 # For temperature calculations done by ParticleTemperaturePropertyUpdater
                 pressure: float = 100.0, thermal_creep: float = 1.0,
                 molar_mass: float = 0.004002602, molar_heat_capacity: float = 12.5,
                 specific_gas_constant: float = 2077.0, specific_heat_capacity: float = 3116.0,
                 speed_of_sound: float = 117.7):
        # Store all the fixed initial properties
        self.dynamic_viscosity = dynamic_viscosity
        self.density = density
        self.temperature = temperature
        self.slip_correction_scale = slip_correction_scale
        self.m_gas = m_gas
        self.pressure = pressure

        # Set the correct slip correction model
        self._init_slip_correction_model(slip_correction_model, temperature)

        # Store additional properties needed for particle temperature calculations
        self.thermal_creep = thermal_creep
        self.molar_mass = molar_mass
        self.molar_heat_capacity = molar_heat_capacity
        self.specific_gas_constant = specific_gas_constant
        self.specific_heat_capacity = specific_heat_capacity
        self.speed_of_sound = speed_of_sound

        # Short-term memoization storage
        self.interpolation_results = {}
        super().__init__(filename)

    def _init_slip_correction_model(self, slip_correction_model: str, temperature: float):
        # Either take the manually set slip correction model,
        if slip_correction_model is not None:
            self.slip_correction_model = slip_correction_model
        else:  # or define it automatically based on the temperature.
            if temperature == 4.0:
                self.slip_correction_model = '4_kelvin'
            elif temperature == 293.15:
                self.slip_correction_model = 'room_temp'
            else:  # Use some heuristic for the model but warn the user that this has been used.
                self.slip_correction_model = '4_kelvin' if 0.0 <= temperature <= 200.0 else 'room_temp'
                print(f"WARNING: Auto-picked slip correction model '{self.slip_correction_model}' based on"
                      f" a temperature here the model might not be applicable. You might need to set the model by hand,"
                      f" or even implement a new model, for optimal results.")

        # In any case, check that the model that has been set is one that is actually implemented
        if self.slip_correction_model not in ['4_kelvin', 'room_temp']:
            raise ValueError(f"Invalid slip correction model: {self.slip_correction_model}")
        if self.slip_correction_model == 'room_temp':
            # Require m_gas for the room_temp model
            if self.m_gas is None:
                raise ValueError("m_gas is a required parameter for the 'room_temp' slip correction model!")

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
        Calculates the drag force using Stokes' law for spherical particles in continuum
        :param relative_velocity: The velocity of the field relative to the particle
        :param pressure: The pressure exerted on the particle at the point
        :param particle: A Particle
        :return: The drag force as a (3,) np.array.
        """
        # Negative pressure or zero pressure is akin to the particle not being in the flow field, force is zero
        # and we can stop considering the particle to be inside this field
        if pressure <= 0:
            return np.zeros(self.number_of_dimensions)

        force_vector = 6 * np.pi * self.dynamic_viscosity * particle.radius * relative_velocity
        return force_vector / self.calc_slip_correction(pressure=pressure, particle=particle)

    def calc_slip_correction(self, pressure: float, particle: SphericalParticle) -> float:
        """
        Calculates the slip correction factor with temperature corrections. Works for models '4_kelvin' at 4K,
        and 'room_temp' at 293.15K.

        The 'room_temp' models is based on the paper:
        D. K. Hutchins, M. H. Harper, R. L. Felder, Slip correction measurements for solid spherical particles
        by modulated dynamic light scattering, Aerosol Sci. Techn. 22 (2) (1995) 202–218. doi:10.1080/02786829408959741.
        """
        if pressure == 0.0:
            return float('inf')  # Correct force to 0 by dividing by infinity. FIXME better way?

        s = 0.0
        if self.slip_correction_model == '4_kelvin':
            # The Sutherland constant for helium is 79.4 at reference temperature of 273.0 K.
            # I took a reference pressure of 1 Pascal, in which the mean free path of helium is 0.01754.
            # see J. Aerosol Sci. 1976 Vol. 7. pp 381-387 by Klaus Willeke
            c_sutherland = 79.4
            ref_temp = 273.0
            knudsen = 0.01754 * 1 / (pressure * particle.radius)\
                      * self.temperature / ref_temp * (1 + c_sutherland / ref_temp)\
                      / (1 + c_sutherland / self.temperature)
            s = 1 + knudsen * (1.246 + (0.42 * np.exp(-0.87 / knudsen)))
        elif self.slip_correction_model == 'room_temp':
            knudsen = self.dynamic_viscosity / (pressure * particle.radius) *\
                      np.sqrt(pi * Boltzmann * self.temperature / (2 * self.m_gas))
            s = 1 + knudsen * (1.231 + 0.4695 * np.exp(-1.1783 / knudsen))
        return s * self.slip_correction_scale

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
