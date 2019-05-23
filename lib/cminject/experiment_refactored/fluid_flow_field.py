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

from cminject.experiment_refactored.base_classes import Field
from cminject.experiment_refactored.basic import SphericalParticle
from scipy.interpolate import RegularGridInterpolator


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
                 molar_mass: float = 0.004002602,
                 # molar_heat_capacity: float = 12.5,
                 # specific_gas_constant: float = 2077.0, specific_heat_capacity: float = 3116.0,
                 # m_gas: float = 6.6e-27, m_gas_mol: float = 0.004002602,
                 # speed_of_sound: float = 117.7,
                 scale_slip: float = 1.0):
        self.density = density
        self.dynamic_viscosity = dynamic_viscosity
        self.temperature = temperature
        self.scale_slip = scale_slip

        self.pressure = pressure  # TODO the pressure of what? averaged?
        self.thermal_creep = thermal_creep
        self.molar_mass = molar_mass

        # Storage to quickly remember and look up particles that are considered to be outside this field
        self.outside_particles = set()

        # Construct the interpolator for the drag force from the HDF5 file passed by file name
        self.filename = filename
        x, y, z, data_grid = self._read_from_hdf5(self.filename)
        self.min_z, self.max_z = np.min(z), np.max(z)
        self.f_drag = RegularGridInterpolator((x, y, z), data_grid)

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.min_z, self.max_z

    def is_particle_inside(self, particle: SphericalParticle) -> bool:
        return self.min_z <= particle.position[2] <= self.max_z\
               and particle.identifier not in self.outside_particles

    def calculate_drag_force(self, particle: SphericalParticle) -> np.array:
        """
        This function calculates the drag force using Stokes' law for spherical particles in continuum
        """
        try:
            f_drag_result = self.f_drag((particle.position[0], particle.position[1], particle.position[2]))
            pressure = f_drag_result[3]
            # Negative pressure or zero pressure is akin to the particle not being in the flow field, force is zero
            # and we can stop considering the particle to be inside this field
            if pressure <= 0:
                self.outside_particles.add(particle.identifier)
                return np.zeros(3)

            relative_velocity = f_drag_result[:3] - particle.velocity
            force_vector = 6 * np.pi * self.dynamic_viscosity * particle.radius * relative_velocity
            return force_vector / self.calculate_slip_correction(pressure=pressure, particle=particle)
        except ValueError:  # Particle is definitely outside the field boundaries, so force is zero
            self.outside_particles.add(particle.identifier)
            return np.zeros(3)

    def calculate_slip_correction(self, pressure: float, particle: SphericalParticle) -> float:
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

    def calculate_acceleration(self, particle: SphericalParticle, time: float) -> np.array:
        return self.calculate_drag_force(particle) / particle.mass

    @staticmethod
    def _read_from_hdf5(filename):
        print(f"Reading in HDF5 file {filename}...")
        from cminject.experiment_refactored import comsol_hdf5_tools
        x, y, z, data_grid = comsol_hdf5_tools.hdf5_to_data_grid(filename)
        print(f"Done reading in HDF5 file {filename}.")
        # Return x/y/z and the data grid. Both are needed to construct a RegularGridInterpolator.
        return x, y, z, data_grid


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""