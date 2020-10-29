#!/usr/bin/env python3
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMInject
#
# Copyright (C) 2018,2020 CFEL Controlled Molecule Imaging group
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any 
# later version. 
#
# If you use this program for scientific work, you should correctly reference it; see the LICENSE.md file for details. 
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Acceleration fields that model flowing fluids.
"""

import logging
import warnings
from abc import ABC
from typing import Tuple

import h5py  # TODO should this code here really be concerned with I/O?
import numpy as np

from cminject.calc import fluid_flow, is_finite
from cminject.particles.spherical import SphericalParticle, ThermallyConductiveSphericalParticle
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

    def get_local_properties(self, phase_space_position: np.array) -> Tuple[float, np.array]:
        """
        Gets the properties of the fluid at the phase-space position (position + velocity) of a particle.

        :param phase_space_position: The phase-space position of the particle.
        :return: A 2-tuple, containing
          * The pressure of the fluid at that point
          * The relative velocity (v_fluid - v_particle)
        """
        data = self.interpolate(phase_space_position[:self.number_of_dimensions])
        relative_velocity = data[:self.number_of_dimensions] - phase_space_position[self.number_of_dimensions:]
        pressure = data[self.number_of_dimensions]
        return pressure, relative_velocity

    def interpolate(self, position: np.array) -> np.array:
        """
        Return the raw results of interpolating on this field at a certain position in space.

        :param position: The position in space to interpolate this field for.
        :return:
          * An np.array containing all interpolated quantities at the given position, or
          * An np.array containing only NaN, if the interpolator raised a ``ValueError``.
        """
        try:
            return self._interpolator(position)
        except ValueError:
            return np.full(self.number_of_dimensions + 1, np.nan)  # vr,vz,p / vx,vy,vz,p / ...

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        return self.interpolate(position)[self.number_of_dimensions] > 0.0

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def _set_cascading(self, attribute_name, passed_arg, h5f, key):
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

    def calculate_acceleration(self, particle: SphericalParticle, time: float) -> np.array:
        """
        Calculates the drag force using Stokes' law for spherical particles in continuum
        """
        pressure, relative_velocity = self.get_local_properties(particle.phase_space_position)

        if pressure > 0 and is_finite(pressure):
            Cc = self.calc_slip_correction(pressure, particle.radius)
            return fluid_flow.a_stokes(relative_velocity, self.dynamic_viscosity, particle.radius, particle.mass, Cc)
        else:
            return self._nan_acceleration

    def calc_slip_correction(self, pressure: float, particle_radius: float) -> float:
        """
        Calculates the slip correction factor with temperature corrections. Works for models '4_kelvin' at 4K,
        and 'room_temp' at 293.15K.

        .. note::
            This method is replaced at runtime with a concrete implementation for the chosen slip correction model,
            to avoid the overhead of choosing at every call.
        """
        raise NotImplementedError(
            "Unknown slip correction model; this method should have been replaced at construction!"
        )

    def _calc_slip_correction_4k(self, pressure: float, particle_radius: float) -> float:
        return fluid_flow.slip_correction_4k(pressure, particle_radius, self.temperature, self.slip_correction_scale)

    def _calc_slip_correction_hutchins(self, pressure: float, particle_radius: float) -> float:
        return fluid_flow.slip_correction_hutchins(
            self.dynamic_viscosity, pressure, particle_radius, self.temperature,
            self.m_gas, self.slip_correction_scale
        )

    def _set_properties_based_on_gas_type(self, gas_type: str, gas_temp: float):
        # TODO this should be some kind of table, maybe along with interpolation, rather than for a very small finite...
        # TODO ...set of temperature/gas combinations
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
        pressure, relative_velocity = self.get_local_properties(particle.phase_space_position)

        # Negative pressure or zero pressure is akin to the particle not being in the flow field, force is zero
        if pressure <= 0:
            return self._zero_acceleration
        elif not np.isfinite(pressure):
            return self._nan_acceleration

        return fluid_flow.a_roth(pressure, relative_velocity, self.m_gas, self.temperature,
                                 particle.temperature, particle.radius, particle.mass)
