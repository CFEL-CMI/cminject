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

import logging
import warnings
from abc import ABC
from typing import Tuple

import h5py
import numba
import numpy as np
from scipy.constants import pi, Boltzmann
from scipy.special import erf

from cminject.definitions.particles.spherical import SphericalParticle
from cminject.definitions.particles.t_conductive_spherical import ThermallyConductiveSphericalParticle

from .regular_grid_interpolation import RegularGridInterpolationField


class DragForceInterpolationField(RegularGridInterpolationField, ABC):
    """
    A base class for any field that calculates a drag force based on a regular interpolation grid field and
    some model. The method that needs to be implemented is `calculate_drag_force`.
    """
    def __init__(self, filename, *args, **kwargs):
        # Short-term memoization storage
        self.interpolation_results = {}
        super().__init__(filename=filename, *args, **kwargs)

    def get_local_properties(self, particle_position: np.array, particle_velocity: np.array) -> Tuple[float, np.array]:
        data = self.interpolate(particle_position)
        relative_velocity = data[:self.number_of_dimensions] - particle_velocity
        pressure = data[self.number_of_dimensions]
        return pressure, relative_velocity

    def interpolate(self, position: np.array) -> np.array:
        try:
            return self._interpolator(position)
        except ValueError:
            return np.zeros(self.number_of_dimensions + 1)  # vr,vz,p / vx,vy,vz,p / ...

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        return self.interpolate(position)[self.number_of_dimensions] > 0.0

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def _set_cascading(self, attribute_name, passed_arg, h5f, key):
        # TODO docstring
        if passed_arg is not None:
            setattr(self, attribute_name, passed_arg)
        elif key in h5f.attrs:
            setattr(self, attribute_name, h5f.attrs[key])
        elif getattr(self, attribute_name, None) is None:
            # couldn't find value and hasn't already been set
            raise ValueError(f"Could not retrieve value for {attribute_name} from HDF5 and function arguments!")


class StokesDragForceField(DragForceInterpolationField):
    """
    A flow field that calculates a drag force exerted on a particle, based on the stokes force and interpolation
    on a grid defined by an HDF5 file like comsol_hdf5_tools.txt_to_hdf5 outputs.
    """
    def __init__(self, filename: str,
                 dynamic_viscosity: float = None, m_gas: float = None, temperature: float = None,
                 slip_correction_model: str = None, slip_correction_scale: float = 1.0,
                 *args, **kwargs):
        # Store all the fixed initial properties
        self.dynamic_viscosity, self.density, self.m_gas, self.temperature, self.slip_correction_scale = [None] * 5
        with h5py.File(filename, 'r') as h5f:
            if 'flow_gas' in h5f.attrs:
                temp = h5f.attrs['flow_temperature']
                self._set_properties_based_on_gas_type(h5f.attrs['flow_gas'], temp)
                if temperature is not None and temperature != temp:
                    warnings.warn(
                        f"WARNING: Given temperature ({temperature}) is not equal to the temperature stored on the "
                        f"field ({temp}). Things like the assumed dynamic viscosity might be off, so "
                        f"you might need to set them manually.")

            # Set all the values that are either stored on the field file or passed explicitly, getting the value
            # to set in a 'cascading' manner, where first the field file is looked at, and the value from there stored,
            # and then the passed argument is looked at, which -- if it is not None -- overrides the value from the
            # field file.
            self._set_cascading('dynamic_viscosity', dynamic_viscosity, h5f, 'flow_dynamic_viscosity')
            self._set_cascading('m_gas', m_gas, h5f, 'flow_gas_mass')
            self._set_cascading('temperature', temperature, h5f, 'flow_temperature')
            self._set_cascading('slip_correction_scale', slip_correction_scale, h5f, 'flow_slip_scale')

        # Set the correct slip correction model
        self.slip_correction_model = None
        self._init_slip_correction_model(slip_correction_model)
        if self.slip_correction_model is None:
            raise ValueError(f"{self.__class__} was unable to set a slip correction model!")

        super().__init__(filename, *args, **kwargs)

    @staticmethod
    @numba.jit('float64[::1](float64, float64[::1], float64, float64, float64, float64, float64[::1])', nopython=True)
    def _a(p, delta_v, mu, r_p, m_p, Cc, azero):
        if p <= 0:
            return azero
        return 6*pi * mu * r_p * delta_v / (m_p * Cc)

    def calculate_acceleration(self, particle: SphericalParticle, time: float) -> np.array:
        """
        Calculates the drag force using Stokes' law for spherical particles in continuum
        """
        pressure, relative_velocity = self.get_local_properties(particle.spatial_position, particle.velocity)
        Cc = self.calc_slip_correction(pressure, particle.radius)
        return self._a(pressure, relative_velocity, self.dynamic_viscosity,
                       particle.radius, particle.mass, Cc, self._zero_acceleration)

    @staticmethod
    @numba.jit(nopython=True)
    def _slip_correction_4k(p, r_p, T, ss):
        # The Sutherland constant for helium is 79.4 at reference temperature of 273.0 K.
        # I took a reference pressure of 1 Pascal, in which the mean free path of helium is 0.01754.
        # see J. Aerosol Sci. 1976 Vol. 7. pp 381-387 by Klaus Willeke
        knudsen = 0.01754 * 1 / (p * r_p) * (T / 273.0) * ((1 + 79.4 / 273.0) / (1 + 79.4 / T))
        s = 1 + knudsen * (1.246 + (0.42 * np.exp(-0.87 / knudsen)))
        return s * ss

    @staticmethod
    @numba.jit(nopython=True)
    def _slip_correction_hutchins(mu, p, r_p, T, m_g, ss):
        # D. K. Hutchins, M. H. Harper, R. L. Felder, Slip correction measurements for solid spherical particles
        # by modulated dynamic light scattering, Aerosol Sci. Techn. 22 (2) (1995) 202–218.
        # doi:10.1080/02786829408959741.
        knudsen = mu / (p * r_p) * np.sqrt(pi * Boltzmann * T / (2 * m_g))
        s = 1 + knudsen * (1.231 + 0.4695 * np.exp(-1.1783 / knudsen))
        return s * ss

    def calc_slip_correction(self, pressure: float, particle_radius: float) -> float:
        """
        Calculates the slip correction factor with temperature corrections. Works for models '4_kelvin' at 4K,
        and 'room_temp' at 293.15K.

        The 'room_temp' models is based on the paper:
        D. K. Hutchins, M. H. Harper, R. L. Felder, Slip correction measurements for solid spherical particles
        by modulated dynamic light scattering, Aerosol Sci. Techn. 22 (2) (1995) 202–218. doi:10.1080/02786829408959741.

        .. todo::
            what model is 4_kelvin based on?

        .. note::
            This method is replaced at runtime with a concrete implementation for the chosen slip correction model,
            to avoid the overhead of choosing at every call.
        """
        raise NotImplementedError("Unknown slip correction model -- "
                                  "this method should have been replaced at construction!")

    def _calc_slip_correction_4k(self, pressure: float, particle_radius: float) -> float:
        if pressure == 0.0:
            return float('inf')  # Correct force to 0 by dividing by infinity. FIXME better way?
        return self._slip_correction_4k(pressure, particle_radius, self.temperature, self.slip_correction_scale)

    def _calc_slip_correction_hutchins(self, pressure: float, particle_radius: float) -> float:
        if pressure == 0.0:
            return float('inf')  # Correct force to 0 by dividing by infinity. FIXME better way?
        return self._slip_correction_hutchins(
            self.dynamic_viscosity, pressure, particle_radius, self.temperature,
            self.m_gas, self.slip_correction_scale
        )

    def _set_properties_based_on_gas_type(self, gas_type: str, gas_temp: float):
        self.temperature = gas_temp

        if gas_type is not None:
            logging.info(f"Automatically setting mu and m_gas for gas type '{gas_type}' at {gas_temp}K.")
        if gas_type == 'He':
            self.m_gas = 6.6e-27
            if gas_temp == 4.0:
                self.dynamic_viscosity = 1.02e-6
            else:
                self.dynamic_viscosity = 1.96e-5
                if gas_temp != 293.15:
                    warnings.warn(
                        f"WARNING: Unknown mu/m for gas type {gas_type} at {gas_temp}. Using "
                        f"approximation for room temperature (293.15K).")
        elif gas_type == 'N':
            self.m_gas = 4.7e-26
            self.dynamic_viscosity = 1.76e-5
            if gas_temp != 293.15:
                warnings.warn(
                    f"WARNING: Unknown mu/m for gas type {gas_type} at {gas_temp}. Using "
                    f"approximation for room temperature (293.15K).")

    def _init_slip_correction_model(self, slip_correction_model: str):
        # Either take the manually set slip correction model,
        if slip_correction_model is not None:
            self.slip_correction_model = slip_correction_model
        else:  # or define it automatically based on the temperature.
            if self.temperature == 4.0:
                self.slip_correction_model = '4_kelvin'
            elif self.temperature == 293.15:
                self.slip_correction_model = 'room_temp'
            else:  # Use some heuristic for the model but warn the user that this has been used.
                self.slip_correction_model = '4_kelvin' if 0.0 <= self.temperature <= 200.0 else 'room_temp'
                warnings.warn(
                    f"WARNING: Auto-picked slip correction model '{self.slip_correction_model}' based on"
                    f" a temperature where the model might not be applicable. You might need to set the model by hand,"
                    f" or even implement a new model, for optimal results.")

        # In any case, check that the model that has been set is one that is actually implemented
        if self.slip_correction_model not in ['4_kelvin', 'room_temp']:
            raise ValueError(f"Invalid slip correction model: {self.slip_correction_model}")
        if self.slip_correction_model == 'room_temp':
            # Require m_gas for the room_temp model
            if self.m_gas is None:
                raise ValueError("m_gas is a required parameter for the 'room_temp' slip correction model!")

        # Swap out the implementation
        if self.slip_correction_model == '4_kelvin':
            self.calc_slip_correction = self._calc_slip_correction_4k
        else:
            self.calc_slip_correction = self._calc_slip_correction_hutchins

        logging.info(f"Set slip correction model to {self.slip_correction_model}.")


class MolecularFlowDragForceField(DragForceInterpolationField):
    """
    A flow field that calculates a drag force exerted on a particle, based on the Epstein force for high velocities and
    interpolation on a grid defined by an HDF5 file like comsol_hdf5_tools.txt_to_hdf5 outputs.
    """
    def __init__(self, filename: str, m_gas: float = None, temperature: float = None, *args, **kwargs):
        self.temperature, self.m_gas = None, None
        with h5py.File(filename, 'r') as h5f:
            self._set_cascading('m_gas', m_gas, h5f, 'flow_gas_mass')
            self._set_cascading('temperature', temperature, h5f, 'flow_temperature')

        super().__init__(filename, *args, **kwargs)

    def calculate_acceleration(self,
                               particle: ThermallyConductiveSphericalParticle,
                               time: float) -> np.array:
        """
        Calculates the drag force using Epstein's law for spherical particles
        in molecular flow with corrections for high velocities
        """
        pressure, relative_velocity = self.get_local_properties(particle.spatial_position, particle.velocity)

        # Negative pressure or zero pressure is akin to the particle not being in the flow field, force is zero
        # and we can stop considering the particle to be inside this field
        if pressure <= 0:
            return np.zeros(self.number_of_dimensions)
        h = self.m_gas / (2*Boltzmann*self.temperature)
        h_ = self.m_gas / (2*Boltzmann*particle.temperature)

        # calculating the force for two different cases
        # 1. Epstein's formula (V<10 m/s)
        # 2. Epstein's formula corrected for high velocities (V>=10 m/s)
        f_spec = np.where(
            abs(relative_velocity) < 10,

            16/3 * pressure * np.sqrt(pi*h) * particle.radius**2 * relative_velocity,

            -pressure * np.sqrt(pi) * particle.radius**2 * (
                -2*np.exp(-h*relative_velocity**2) * np.sqrt(h) * relative_velocity * (1+2*h*relative_velocity**2) +
                np.sqrt(pi) * (1 - 4 * h * relative_velocity**2 - 4 * h**2 * relative_velocity**4) *
                erf(np.sqrt(h) * relative_velocity)
            ) / (2 * h * relative_velocity**2)
        )
        force_vector = f_spec + 1.8 / 3 * pressure * pi**(3/2) * h/np.sqrt(h_) * particle.radius**2 * relative_velocity
        acceleration = force_vector / particle.mass
        return acceleration

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
