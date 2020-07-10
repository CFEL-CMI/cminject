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

from typing import Callable, Tuple

import numpy as np

from cminject.definitions.particles.base import Particle
from cminject.definitions.util import infinite_interval

from .base import Field


class FunctionField(Field):
    """
    A simple stateless Field based on a function to calculate the acceleration.

    Implies infinite extent of the field, i.e. a particle is never considered to be outside the field's boundary.
    The number of dimensions of this field is not validated, meaning the given function must be able to handle
    particles of the experiment dimensionality.
    """
    def calculate_acceleration(self, particle: Particle, time: float) -> np.array:
        return self.function(particle, time)

    def __init__(self, function: Callable[[Particle, float], np.array]):
        self.function = function

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
