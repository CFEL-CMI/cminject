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

from typing import List, Type, Union, Dict, Any

import numpy as np

from cminject.base import Source
from cminject.particles.spherical import SphericalParticle
from cminject.utils.global_config import GlobalConfig, ConfigSubscriber, ConfigKey

Distribution = Union[float, Dict[str, float]]


class VariableDistributionSource(Source, ConfigSubscriber):
    """
    A Source for particles that allows variable distributions for each dimension of the initial phase space.
    A "Distribution" is a dictionary describing a specific kind of distribution, with a 'kind' key that must match
    one of the following distributions, and further keys describing parameters for the given distribution:

    - gaussian: A gaussian distribution with keys 'mu' and 'sigma'. (like np.random.normal)
    - linear: A deterministically linear distribution with keys 'min' and 'max' (like np.linspace)
    - uniform: A random uniform distribution with keys 'min' and 'max' (like np.random.uniform)
    - radial_gaussian: A 1D-projected 2D gaussian distribution. Useful for a radial dimension.
    - radial_linear: An approximation of a 1D-projected 2D deterministically linear distribution.
        Useful for a radial dimension.
    - radial_uniform: A 1D-projected 2D random uniform distribution. Useful for a radial dimension.
    """
    def __init__(self,
                 number_of_particles: int,
                 position: List[Distribution],
                 velocity: List[Distribution],
                 radius: Distribution,
                 density: float,
                 seed=None,
                 subclass: Type[SphericalParticle] = SphericalParticle,
                 subclass_kwargs: Dict[Any, Any] = None):
        self.number_of_particles = number_of_particles
        self.position = position
        self.velocity = velocity
        self.radius = radius
        self.rho = density
        self.seed = seed
        self.subclass = subclass
        self.subclass_kwargs = subclass_kwargs or {}
        self.number_of_dimensions = len(self.position)

        # Check that position and velocity descriptions are equal in length (dimensionality)
        if len(self.position) != len(self.velocity):
            raise ValueError(
                f"Position and velocity descriptions must match in size! "
                f"Given position description was {len(self.position)}-dimensional, "
                f"velocity description was {len(self.velocity)}-dimensional."
            )
        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            if value != self.number_of_dimensions:
                raise ValueError(
                    f"Incompatible number of dimensions for {self}: {value}, "
                    f"the position description passed at construction was {self.number_of_dimensions}-dimensional."
                )

    def _generate(self, dist: Distribution):
        if type(dist) is float:
            return np.repeat(dist, self.number_of_particles)

        if dist['kind'] == 'gaussian':
            mu, sigma = dist['mu'], dist['sigma']
            return np.random.normal(mu, sigma, self.number_of_particles)
        elif dist['kind'] == 'linear':
            l, r = dist['min'], dist['max']
            return np.linspace(l, r, self.number_of_particles)
        elif dist['kind'] == 'radial_gaussian':
            mu, sigma = dist['mu'], dist['sigma']
            random_x = np.random.normal(mu, sigma, self.number_of_particles)
            random_y = np.random.normal(0, sigma, self.number_of_particles)
            return np.sqrt(random_x**2 + random_y**2)
        elif dist['kind'] == 'radial_linear':
            l, r = dist['min'], dist['max']
            lin = np.linspace(l, r, self.number_of_particles)
            h = 1/np.sqrt(np.abs(lin[np.where(lin != 0)]))
            c = np.cumsum(h)
            return l + (c*(r-l) / np.sum(h))
        elif dist['kind'] == 'uniform':
            l, r = dist['min'], dist['max']
            return np.random.uniform(l, r, self.number_of_particles)
        elif dist['kind'] == 'radial_uniform':
            l, r = dist['min'], dist['max']
            random_x = np.random.uniform(l, r, self.number_of_particles)
            random_y = np.random.uniform(l, r, self.number_of_particles)
            return np.sqrt(random_x**2 + random_y**2)
        else:
            raise ValueError(f"Unknown or unspecified kind in distribution specification: {dist}")

    def generate_particles(self, start_time: float = 0.0):
        position = np.array([self._generate(pdist) for pdist in self.position]).transpose()
        velocity = np.array([self._generate(vdist) for vdist in self.velocity]).transpose()
        r = self._generate(self.radius)

        particles = []
        for i in range(self.number_of_particles):
            inst = self.subclass(
                identifier=i,
                position=position[i],
                velocity=velocity[i],
                start_time=start_time,
                radius=r[i],
                rho=self.rho,
                **self.subclass_kwargs
            )
            particles.append(inst)
        return particles

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
