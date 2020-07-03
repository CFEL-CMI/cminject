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

from abc import ABC
from typing import Tuple, Any

import numpy as np

from cminject.tools.structured_txt_hdf5_tools import hdf5_to_data_grid
from cminject.utils.interpolation import get_regular_grid_interpolator
from cminject.global_config import GlobalConfig, ConfigSubscriber, ConfigKey
from .base import Field


class RegularGridInterpolationField(Field, ConfigSubscriber, ABC):
    """
    A generic base class for fields that calculate a force based on interpolation within some n-dimensional
    regular grid. This class provides an `interpolate` method; the actual calculation of acceleration has to be
    defined in a subclass and should use the result of this method.
    """
    def __init__(self, filename: str, offset: np.array = None, *args, **kwargs):
        """
        The constructor for RegularGridInterpolationField.

        :param filename: The filename of an HDF5 file in the format as constructed by `tools/txt_to_hdf5`.
        :param offset: The positional offset of the field as a (d,)-ndarray, where d is the number of spatial dimensions

        .. todo:: Provide a "formal" definition of the Field-file format in documentaiton, not only through an implementation
           And then refer to this file format by its name, not by its "implemented somewhere" specification...

        """
        # Construct the interpolator from the HDF5 file passed by file name
        self.filename = filename
        data_index, data_grid = hdf5_to_data_grid(self.filename)
        if offset is not None:
            for i, off in enumerate(offset):
                data_index[i] += off

        self.number_of_dimensions = len(data_index)
        self._zero_acceleration = np.zeros(self.number_of_dimensions)
        self._interpolator = get_regular_grid_interpolator(tuple(data_index), data_grid)

        # Assuming that Z is always the last dimension, construct the Z boundary from it
        minima = list(np.min(d) for d in data_index)
        maxima = list(np.max(d) for d in data_index)
        self._z_boundary = (minima[self.number_of_dimensions - 1], maxima[self.number_of_dimensions - 1])
        self.grid_boundary = np.array(list(zip(minima, maxima)))
        super().__init__(*args, **kwargs)

        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            if value != self.number_of_dimensions:
                raise ValueError(
                    f"This field was constructed from a {self.number_of_dimensions}-dimensional grid and is "
                    f"thus incompatible with a {value}-dimensional simulation!"
                )

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        return np.all((position >= self.grid_boundary[:, 0]) * (position <= self.grid_boundary[:, 1]))

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
