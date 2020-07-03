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

from typing import List, Tuple, Any

import numpy as np

from .simple_z import SimpleZBoundary
from cminject.global_config import ConfigSubscriber, GlobalConfig, ConfigKey


class CuboidBoundary(SimpleZBoundary, ConfigSubscriber):
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
        return np.all(
            (self.intervals[:, 0] <= position) * (position <= self.intervals[:, 1])
        )

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
