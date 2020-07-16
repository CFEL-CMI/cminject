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

from functools import cached_property

import numpy as np

from cminject.base import Particle


class SphericalParticle(Particle):
    """
    A simple spherical particle that has a radius and a density (rho).
    """
    def __init__(self, *args, radius: float, rho: float, **kwargs):
        self.radius = radius
        self.rho = rho
        super().__init__(*args, **kwargs)

    @property
    def constant_properties(self):
        return super().constant_properties + [
            ('radius', np.float64),
            ('rho', np.float64)
        ]

    @property
    def tracked_properties(self):
        return super().tracked_properties

    @cached_property
    def mass(self):
        return self.rho * 4 / 3 * np.pi * (self.radius ** 3)


class ThermallyConductiveSphericalParticle(SphericalParticle):
    """
    Lacking a better name (so far), this class extends on the SphericalParticle with a thermal conductivity
    that is required by e.g. the photophoretic force and for thermal calculations.
    """
    def __init__(self, *args, temperature: float, thermal_conductivity: float = None, specific_heat: float = None,
                 **kwargs):
        self.thermal_conductivity = thermal_conductivity
        self.temperature = temperature
        self.specific_heat = specific_heat
        super().__init__(*args, **kwargs)

    @property
    def constant_properties(self):
        return super().constant_properties + [
            ('thermal_conductivity', np.float64),
            ('specific_heat', np.float64)
        ]

    @property
    def tracked_properties(self):
        return super().tracked_properties + [
            ('temperature', np.float64)
        ]


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End: