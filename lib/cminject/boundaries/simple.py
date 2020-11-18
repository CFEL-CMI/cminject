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

"""A collection of simple implementations of :class:`cminject.base.Boundary`."""

from typing import Tuple, List, Any

import numpy as np

from cminject.base import Boundary
from cminject.utils import infinite_interval
from cminject.utils.global_config import ConfigSubscriber, GlobalConfig, ConfigKey


class SimpleZBoundary(Boundary):
    """
    A boundary defined in terms of an interval in the Z dimension (the last spatial dimension).
    This type of boundary should be considered infinite with respect to all other spatial dimensions.
    """
    def __init__(self, z_minmax: Tuple[float, float]):
        z_min, z_max = z_minmax
        if z_min > z_max:
            raise ValueError("z_min must be < z_max!")
        self.z_min = z_min
        self.z_max = z_max
        self._z_boundary = (z_min, z_max)
        super().__init__()

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, position: np.array, time: float):
        """
        True if the particle's position is inside the defined interval in the Z dimension (the last spatial dimension).
        All other coordinates are ignored, so this boundary should be considered infinite with respect to those.

        :param position: The particle's position (ignored apart from the last component
        :param time: The current time (ignored)
        """
        return self.z_min <= position[-1] <= self.z_max


class InfiniteBoundary(Boundary):
    """
    An infinite boundary, enclosing all of space.
    """
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def is_particle_inside(self, position: float, time: float) -> bool:
        return True


class CuboidBoundary(SimpleZBoundary, ConfigSubscriber):
    """
    An n-dimensional cuboid (line/rectangle/3D cuboid...) that encloses a region of space.
    """
    def __init__(self, intervals: List[Tuple[float, float]]):
        self.intervals = np.array(intervals)
        super().__init__(intervals[-1])  # last dimension is assumed to be Z

        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            number_of_dimensions, interval_dimensions = value, len(self.intervals)
            if number_of_dimensions != interval_dimensions:
                raise ValueError(
                    f"Incompatible number of dimensions: {number_of_dimensions}, "
                    f"the cuboid is defined within {interval_dimensions}-dimensional space!"
                )

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        """
        True if the particle's position is inside the defined cuboid region.

        :param position: The particle's position.
        :param time: The current time (ignored)
        """
        return np.all(
            (self.intervals[:, 0] <= position) * (position <= self.intervals[:, 1])
        )
