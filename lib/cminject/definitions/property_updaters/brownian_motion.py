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

import numba
import numpy as np
from cminject.definitions.particles.spherical import SphericalParticle
from cminject.definitions.particles.t_conductive_spherical import ThermallyConductiveSphericalParticle
from scipy.constants import Boltzmann, pi

from cminject.definitions.fields.fluid_flow import StokesDragForceField, MolecularFlowDragForceField
from cminject.global_config import GlobalConfig, ConfigSubscriber, ConfigKey
from .base import PropertyUpdater


class BrownianMotionPropertyUpdater(PropertyUpdater, ConfigSubscriber):
    """
    Models brownian motion for a Stokes drag force field based on the paper:
    A. Li, G. Ahmadi, Dispersion and deposition of spherical particles from point sources in a turbulent channel flow,
    Aerosol Sci. Techn. 16 (24) (1992) 209–226. doi:10.1080/02786829208959550

    .. note::
        This property updater assumes that when the number of spatial dimensions is 2, the problem is radially
        symmetrical. If you need an x/y 2D simulation rather than a r/z simulation, derive a new class.
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

    def update(self, particle: SphericalParticle, time: float) -> bool:
        data = self.field.interpolate(particle.spatial_position)
        pressure = data[self.number_of_dimensions]
        if pressure >= 0.0:
            cunningham = self.field.calc_slip_correction(pressure, particle.radius)
            s0 = 216 * self.field.dynamic_viscosity * Boltzmann * self.field.temperature /\
                 (pi**2 * (2 * particle.radius)**5 * particle.rho**2 * cunningham)
            a = np.random.normal(0.0, 1.0, self.number_of_dimensions) * np.sqrt(pi * s0 / self.dt)

            position = particle.spatial_position + (0.5 * a * self.dt**2)
            velocity = particle.velocity + (a * self.dt)
            particle.position = np.concatenate((position, velocity))
            return True
        return False


class BrownianMotionMolecularFlowPropertyUpdater(PropertyUpdater, ConfigSubscriber):
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

    def update(self, particle: ThermallyConductiveSphericalParticle, time: float) -> bool:
        pressure = self.field.interpolate(particle.spatial_position)[self.number_of_dimensions]
        if pressure >= 0.0:
            h = self.field.m_gas / (2 * Boltzmann * self.field.temperature)
            h_ = self.field.m_gas / (2 * Boltzmann * particle.temperature)

            s0 = (16/3 + 2/3*pi*np.sqrt(h/h_)) * np.sqrt(pi/h) * pressure * self.field.m_gas * particle.radius**2
            if self.number_of_dimensions == 2:  # TODO this is just for radial symmetry
                a = np.random.normal(0.0, 1.0, 3) * np.sqrt(s0 / self.dt) / particle.mass

                new_x = particle.position[0] + (0.5 * a[0] * self.dt**2)
                new_vx = particle.position[2] + (a[0] * self.dt)
                r_sign = np.sign(new_x)
                vr_sign = np.sign(new_vx)

                particle.position[0] = r_sign * np.sqrt(new_x**2 + (0.5 * a[1] * self.dt**2)**2)
                particle.position[1] = particle.position[1] + (0.5 * a[2] * self.dt**2)
                particle.position[2] = vr_sign * np.sqrt(new_vx**2 + (a[1] * self.dt)**2)
                particle.position[3] = particle.position[3] + (a[2] * self.dt)
            else:
                a = np.random.normal(0.0, 1.0, self.number_of_dimensions) * np.sqrt(s0 / self.dt) / particle.mass
                position = particle.spatial_position + (0.5 * a * self.dt**2)
                velocity = particle.velocity + (a * self.dt)
                particle.position = np.concatenate([position, velocity])

        return True

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
