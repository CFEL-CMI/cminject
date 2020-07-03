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

from typing import List

import numpy as np

from .base import Particle


class SphericalParticle(Particle):
    """
    A simple spherical particle that has a radius and a density (rho).
    """
    def __init__(self, *args, radius: float, rho: float, **kwargs):
        super().__init__(*args, **kwargs)
        self.radius = radius
        self.rho = rho
        self.mass = self.rho * 4 / 3 * np.pi * (self.radius ** 3)

    @property
    def properties(self) -> np.array:
        return np.array([
            self.identifier,
            self.time_of_flight,
            self.radius,
            self.mass
        ])

    @property
    def properties_description(self) -> List[str]:
        return ['ID', 't', 'r', 'm']

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
