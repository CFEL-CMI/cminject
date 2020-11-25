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

"""Actions that affect a Particle's ``temperature`` in some way."""

from typing import Any

import numpy as np
from scipy.constants import Boltzmann

from cminject.base import Action
from cminject.fields.fluid_flow import MolecularFlowDragForceField
from cminject.particles.spherical import ThermallyConductiveSphericalParticle
from cminject.utils.global_config import GlobalConfig, ConfigSubscriber, ConfigKey


class MolecularFlowUpdateTemperature(Action, ConfigSubscriber):
    """
    Updates the temperature of a ThermallyConductiveSphericalParticle, based on the pressure exerted by a
    MolecularFlowDragForceField at the point the particle is at, the particle's specific heat and mass,
    and the chosen time step. The formulae are described in Nils Roth's paper, https://arxiv.org/abs/2006.10652
    """
    def __init__(self, field: MolecularFlowDragForceField):
        self.field = field
        self.dt = None
        self.number_of_dimensions = None
        GlobalConfig().subscribe(self, [ConfigKey.NUMBER_OF_DIMENSIONS, ConfigKey.TIME_STEP])

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            self.number_of_dimensions = value
        if key is ConfigKey.TIME_STEP:
            self.dt = value

    def __call__(self, particle: ThermallyConductiveSphericalParticle, time: float) -> bool:
        pressure = self.field.interpolate(particle.position)[self.number_of_dimensions]
        h = self.field.m_gas / (2 * Boltzmann * self.field.temperature)
        h_ = self.field.m_gas / (2 * Boltzmann * particle.temperature)
        deltaE = 4 * pressure * np.sqrt(np.pi) * particle.radius**2 * (np.sqrt(h)/h_ - 1/np.sqrt(h))
        deltaT = deltaE / (particle.specific_heat * particle.mass)

        particle.temperature -= deltaT * self.dt
        return False
