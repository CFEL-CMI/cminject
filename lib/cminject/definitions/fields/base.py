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

from abc import ABC, abstractmethod

import numpy as np

from cminject.definitions.base import NDimensional, ZBounded
from cminject.definitions.particles.base import Particle


class Field(NDimensional, ZBounded, ABC):
    """
    A Field interacting with Particle objects.
    """
    @abstractmethod
    def calculate_acceleration(self,
                               particle: Particle,
                               time: float) -> np.array:
        """
        Calculates an acceleration for one particle based on the particle's current properties and the current time.
        This acceleration will be integrated for in each time step and thus "applied" to the particle.

        :param particle: The Particle to calculate this Field's acceleration for.
        :param time: The time to calculate the acceleration for.
        :return: A (n,)-shaped numpy array describing the acceleration exerted on the particle. n is the number of
            spatial dimensions of the experiment.
        """
        pass

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
