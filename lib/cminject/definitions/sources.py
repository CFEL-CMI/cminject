#!/usr/bin/env python
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

from typing import List, Type, Dict

import numpy as np
from cminject.definitions.base import Source
from cminject.definitions.particles import SphericalParticle


Distribution = Dict  # An appropriate type for variable distribution descriptions


class VariableDistributionSource(Source):
    def __init__(self,
                 number_of_particles: int,
                 position: List[Distribution],
                 velocity: List[Distribution],
                 radius: Distribution,
                 rho: float,
                 seed=None,
                 randomly_rotate_around_z: bool = False,
                 subclass: Type[SphericalParticle] = SphericalParticle,
                 **subclass_kwargs):
        self.number_of_particles = number_of_particles
        self.position = position
        self.velocity = velocity
        self.radius = radius
        self.rho = rho
        self.seed = seed
        self.randomly_rotate_around_z = randomly_rotate_around_z
        self.subclass = subclass
        self.subclass_kwargs = subclass_kwargs
        self.number_of_dimensions = len(self.position)

        # Check validity of distribution kinds
        for dist in self.position + self.velocity + [self.radius]:
            if 'kind' not in dist or dist['kind'] not in ['gaussian', 'linear']:
                raise ValueError(f"Unknown or unspecified kind in distribution specification: {dist}")
        # Check that position and velocity descriptions are equal in length (dimensionality)
        if len(self.position) != len(self.velocity):
            raise ValueError(
                f"Position and velocity descriptions must match in size! "
                f"Given position description was {len(self.position)}-dimensional, "
                f"velocity description was {len(self.velocity)}-dimensional."
            )

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"Incompatible number of dimensions: {number_of_dimensions}, "
                f"the position description passed at construction was {self.number_of_dimensions}-dimensional."
            )

    def _generate(self, dist: Distribution):
        if dist['kind'] == 'gaussian':
            mu, sigma = dist['mu'], dist['sigma']
            return np.random.normal(mu, sigma, self.number_of_particles)
        elif dist['kind'] == 'linear':
            l, r = dist['min'], dist['max']
            return np.linspace(l, r, self.number_of_particles)
        elif dist['kind'] == 'radial_gaussian':
            mu, sigma = dist['mu'], dist['sigma']
            random_x=np.random.normal(mu, sigma, self.number_of_particles)
            random_y=np.random.normal(0, sigma, self.number_of_particles)
            return np.sqrt(random_x**2 + random_y**2)
        elif dist['kind'] == 'radial_linear':
            l, r = dist['min'], dist['max']
            lin=np.linspace(l,r,self.number_of_particles)
            h=1/np.abs(lin[np.where(lin!=0)])**0.5
            c=np.cumsum(h)
            return l+(c*(r-l)/np.sum(h))

    def _rotate_around_z(self, positions):
        size = positions.shape[0]
        rotated = (positions[:, 0] + positions[:, 1] * 1j) * np.exp(np.random.random(size) * 2 * np.pi * 1j)
        return np.array([rotated.real, rotated.imag, positions[:, 2]]).transpose()

    def generate_particles(self, start_time: float = 0.0):
        position = np.array([self._generate(pdist) for pdist in self.position]).transpose()
        if self.randomly_rotate_around_z:
            position = self._rotate_around_z(position)

        velocity = np.array([self._generate(vdist) for vdist in self.velocity]).transpose()
        r = self._generate(self.radius)

        particles = []
        for i in range(self.number_of_particles):
            inst = self.subclass(
                identifier=i,
                position=np.concatenate([position[i], velocity[i]]),
                start_time=start_time,
                radius=r[i],
                rho=self.rho,
                **self.subclass_kwargs
            )
            particles.append(inst)
        return particles


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""