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

from typing import List, Type, Dict, Any

import numpy as np

from cminject.base import Source
from cminject.particles.spherical import SphericalParticle
from cminject.utils.distributions import Distribution, constant
from cminject.utils.global_config import GlobalConfig, ConfigSubscriber, ConfigKey


class VariableDistributionSource(Source, ConfigSubscriber):
    """
    A Source for particles that allows variable distributions for each dimension of the initial phase space.
    A distribution is a subclass of :class:`cminject.utils.distributions.Distribution`.
    """
    def __init__(self,
                 number_of_particles: int,
                 position: List[Distribution],
                 velocity: List[Distribution],
                 radius: Distribution,
                 density: Distribution,
                 seed=None,
                 particle_class: Type[SphericalParticle] = SphericalParticle,
                 particle_kwargs: Dict[Any, Any] = None):
        self.number_of_particles = number_of_particles
        self.position = position
        self.velocity = velocity
        self.radius = radius
        self.density = density
        self.seed = seed
        self.particle_class = particle_class
        self.particle_kwargs = particle_kwargs or {}
        self.number_of_dimensions = len(self.position)

        # Check that position and velocity descriptions are equal in length (dimensionality)
        m, n = len(self.position), len(self.velocity)
        assert m == n, f"Mismatching position ({m}) / velocity ({n}) dimensionalities!"
        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            if value != self.number_of_dimensions:
                raise ValueError(
                    f"Incompatible number of dimensions for {self}: {value}, "
                    f"the position description passed at construction was {self.number_of_dimensions}-dimensional."
                )

    def generate_particles(self, start_time: float = 0.0):
        n = self.number_of_particles
        position = np.array([pdist.generate(n) for pdist in self.position]).transpose()
        velocity = np.array([vdist.generate(n) for vdist in self.velocity]).transpose()
        r = self.radius.generate(n)
        rho = self.density.generate(n)

        particles = [
            self.particle_class(
                identifier=i,
                position=position[i],
                velocity=velocity[i],
                start_time=start_time,
                radius=r[i],
                rho=rho[i],
                **self.particle_kwargs
            )
            for i in range(n)
        ]
        return particles
