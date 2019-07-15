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

from abc import ABC
from typing import Tuple

import numpy as np
from cminject.definitions.base import Field, Particle
from cminject.tools.structured_txt_hdf5_tools import hdf5_to_data_grid
from scipy.interpolate import RegularGridInterpolator


class RegularGridInterpolationField(Field, ABC):
    """
    A generic base class for fields that calculate a force based on interpolation within some n-dimensional
    regular grid. This class provides an `interpolate` method; the actual calculation of acceleration has to be
    defined in a subclass and should use the result of this method.
    """
    def __init__(self, filename: str):
        """
        The constructor for RegularGridInterpolationField.

        :param filename: The filename of an HDF5 file in the format as constructed by `tools/txt_to_hdf5`.
        """
        # Construct the interpolator from the HDF5 file passed by file name
        self.filename = filename
        data_index, data_grid = hdf5_to_data_grid(self.filename)
        self.number_of_dimensions = len(data_index)
        self.interpolator = RegularGridInterpolator(tuple(data_index), data_grid)

        # Assuming that Z is always the last dimension, construct the Z boundary from it
        minima = list(np.min(d) for d in data_index)
        maxima = list(np.max(d) for d in data_index)
        self._z_boundary = (minima[self.number_of_dimensions - 1], maxima[self.number_of_dimensions - 1])
        self.grid_boundary = np.array(list(zip(minima, maxima)))
        super().__init__()

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, particle: Particle, time: float) -> bool:
        return np.all(particle.spatial_position >= self.grid_boundary[:, 0]) and\
               np.all(particle.spatial_position <= self.grid_boundary[:, 1])

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"This field was constructed from a {self.number_of_dimensions}-dimensional grid and is "
                f"thus incompatible with a {number_of_dimensions}-dimensional simulation!"
            )


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
