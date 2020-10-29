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

"""
Actions that model Brownian motion by applying a random force (with a reasonable spectral intensity) after every
integration step.
"""

from typing import Any

from cminject.base import Action
from cminject.calc import fluid_flow as calc
from cminject.fields.fluid_flow import StokesDragForceField, MolecularFlowDragForceField
from cminject.particles.spherical import SphericalParticle, ThermallyConductiveSphericalParticle
from cminject.utils.global_config import GlobalConfig, ConfigSubscriber, ConfigKey


class StokesBrownianMotionStep(Action, ConfigSubscriber):
    """
    Models brownian motion for a Stokes drag force field based on the paper:
    A. Li, G. Ahmadi, Dispersion and deposition of spherical particles from point sources in a turbulent channel flow,
    Aerosol Sci. Techn. 16 (24) (1992) 209–226. doi:10.1080/02786829208959550
    """
    def __init__(self, field: StokesDragForceField):
        self.field = field
        self.number_of_dimensions = None
        self.dt = None
        GlobalConfig().subscribe(self, [ConfigKey.NUMBER_OF_DIMENSIONS, ConfigKey.TIME_STEP])

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            self.number_of_dimensions = value
        elif key is ConfigKey.TIME_STEP:
            self.dt = value

    def __call__(self, particle: SphericalParticle, time: float) -> bool:
        pressure = self.field.interpolate(particle.position)[self.number_of_dimensions]
        if pressure >= 0.0:
            cunningham = self.field.calc_slip_correction(pressure, particle.radius)
            a = calc.a_brown_stokes(
                self.field.dynamic_viscosity, self.field.temperature,
                particle.radius, particle.rho, cunningham, self.dt, self.number_of_dimensions
            )
            particle.position += 0.5 * a * self.dt**2
            particle.velocity += a * self.dt
            return True
        return False


class MolecularFlowBrownianMotionStep(Action, ConfigSubscriber):
    """
    Models brownian motion based on the fluctuation-dissipation-theorem
    and a numerical representation of the delta function
    """
    def __init__(self, field: MolecularFlowDragForceField):
        self.field = field
        self.number_of_dimensions = 2
        self.dt = None
        GlobalConfig().subscribe(self, [ConfigKey.NUMBER_OF_DIMENSIONS, ConfigKey.TIME_STEP])

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            if value != self.number_of_dimensions:
                raise ValueError(f"{self} requires a {self.number_of_dimensions}-dimensional simulation")
        elif key is ConfigKey.TIME_STEP:
            self.dt = value

    def __call__(self, particle: ThermallyConductiveSphericalParticle, time: float) -> bool:
        pressure = self.field.interpolate(particle.position)[self.number_of_dimensions]
        if pressure >= 0.0:
            a = calc.a_brown_roth(
                pressure, self.field.m_gas, self.field.temperature,
                particle.temperature, particle.radius, particle.mass, self.dt, self.number_of_dimensions
            )
            particle.position += 0.5 * a * self.dt**2
            particle.velocity += a * self.dt
            return True
        return False
