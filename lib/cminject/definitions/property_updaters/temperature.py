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
from typing import Any

import numpy as np
from scipy.constants import Boltzmann

from cminject.definitions.property_updaters.base import PropertyUpdater
from cminject.definitions.fields.fluid_flow import MolecularFlowDragForceField
from cminject.definitions.particles.t_conductive_spherical import ThermallyConductiveSphericalParticle
from cminject.global_config import GlobalConfig, ConfigSubscriber, ConfigKey


class ParticleTemperaturePropertyUpdater(PropertyUpdater, ConfigSubscriber):
    def __init__(self, field: MolecularFlowDragForceField, dt: float):
        self.field = field
        self.dt = dt
        self.number_of_dimensions = None
        GlobalConfig().subscribe(self, [ConfigKey.NUMBER_OF_DIMENSIONS])

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            self.number_of_dimensions = value

    def update(self, particle: ThermallyConductiveSphericalParticle, time: float) -> bool:
        pressure = self.field.interpolate(particle.position)[self.number_of_dimensions]
        h = self.field.m_gas / (2 * Boltzmann * self.field.temperature)
        h_ = self.field.m_gas / (2 * Boltzmann * particle.temperature)
        deltaE = 4 * pressure * np.sqrt(np.pi) * particle.radius**2 * (np.sqrt(h)/h_ - 1/np.sqrt(h))
        deltaT = deltaE / (particle.specific_heat * particle.mass)

        particle.temperature += deltaT * self.dt
        return False

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
