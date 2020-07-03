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

import numpy as np

from cminject.definitions.particles.base import Particle
from .base import PropertyUpdater


class TrajectoryPropertyUpdater(PropertyUpdater):
    """
    A simple property updater to append to the trajectory of the particle after each time step.
    """
    def update(self, particle: Particle, time: float) -> bool:
        particle.trajectory.append(
            np.append(particle.position, time)
        )
        return False

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
