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

from typing import Tuple, List

import numpy as np
from cminject.definitions.base import Boundary, infinite_interval
from cminject.definitions.fields.regular_grid_interpolation_field import RegularGridInterpolationField
from cminject.definitions.fields.stark_field import B_StarkField


class SimpleZBoundary(Boundary):
    def __init__(self, z_minmax: Tuple[float, float]):
        z_min, z_max = z_minmax
        if z_min > z_max:
            raise ValueError("z_min must be < z_max!")
        self.z_min = z_min
        self.z_max = z_max
        self._z_boundary = (z_min, z_max)

        self.number_of_dimensions = 3
        super().__init__()

    def set_number_of_dimensions(self, number_of_dimensions: int):
        self.number_of_dimensions = number_of_dimensions

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, position: np.array, time: float):
        return self.z_min <= position <= self.z_max


class CuboidBoundary(SimpleZBoundary):
    def __init__(self, intervals: List[Tuple[float, float]]):
        self.intervals = np.array(intervals)
        super().__init__(intervals[-1])  # last dimension is assumed to be Z

    def set_number_of_dimensions(self, number_of_dimensions: int):
        interval_dimensions = len(self.intervals)
        if number_of_dimensions != interval_dimensions:
            raise ValueError(
                f"Incompatible number of dimensions: {number_of_dimensions}, "
                f"the cuboid is defined within {interval_dimensions}-dimensional space"
            )
        self.number_of_dimensions = number_of_dimensions

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        return np.all(
            (self.intervals[:, 0] <= position) * (position <= self.intervals[:, 1])
        )


class InfiniteBoundary(Boundary):
    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def is_particle_inside(self, position: float, time: float) -> bool:
        return True


class GridFieldBasedBoundary(Boundary):
    def set_number_of_dimensions(self, number_of_dimensions: int):
        d = number_of_dimensions
        fd = self.field.number_of_dimensions
        if d != fd:
            raise ValueError(f"Dimensionality {d} does not match associated field's dimensionality {fd}!")

    def __init__(self, field: RegularGridInterpolationField):
        self.field = field

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.field.z_boundary

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        return self.field.is_particle_inside(position, time)

class StarkBoundary(Boundary):
    """
    We need to construct a Stark Boundary because the field files expected from Comsol are two dimensional, while the
    experiment is done in 3 dimensions. The z-dimensions are obtained from the user.
    """
    # still need to edit this to take into account the rods of the field
    # Probably I can contruct this class by inheriting GridFieldBasedBoundary and SimpleZBoundary
    def set_number_of_dimensions(self, number_of_dimensions: int):
        d = number_of_dimensions

    def __init__(self, field: B_StarkField, z_minmax: Tuple[float, float]):
        z_min, z_max = z_minmax
        if z_min > z_max:
            raise ValueError("z_min must be < z_max!")
        self.z_min = z_min
        self.z_max = z_max
        self.field = field

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, position: np.array, time: float) -> bool:

        return  self.field.is_particle_inside(position, time) & self.z_min <= position <= self.z_max

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
