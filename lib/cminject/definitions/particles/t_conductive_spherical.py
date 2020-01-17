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

from .spherical import SphericalParticle


class ThermallyConductiveSphericalParticle(SphericalParticle):
    """
    Lacking a better name (so far), this class extends on the SphericalParticle with a thermal conductivity
    that is required by e.g. the photophoretic force and for thermal calculations.
    """
    def __init__(self, *args,
                 thermal_conductivity: float, temperature: float, specific_heat: float,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.thermal_conductivity = thermal_conductivity
        self.temperature = temperature
        self.specific_heat = specific_heat

    @property
    def properties(self) -> np.array:
        return np.concatenate([super().properties, [self.temperature]])

    @property
    def properties_description(self) -> List[str]:
        return super().properties_description + ['T']

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
