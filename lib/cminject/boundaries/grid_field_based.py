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

from typing import Tuple

import numpy as np

from cminject.base import Boundary
from cminject.fields.regular_grid_interpolation import RegularGridInterpolationField


class GridFieldBasedBoundary(Boundary):
    def __init__(self, field: RegularGridInterpolationField):
        self.field = field

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.field.z_boundary

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        return self.field.is_particle_inside(position, time)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
