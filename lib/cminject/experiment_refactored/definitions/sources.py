"""
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
"""

from typing import Tuple, List, Type

import numpy as np
from cminject.experiment_refactored.definitions.base import Source
from cminject.experiment_refactored.definitions.particles import SphericalParticle


class GaussianDistributedSphericalSource(Source):
    """
    A Source generating a Gaussian distribution of SphericalParticles wrt position, velocity and radius.

    By passing `subclass` and any additional arguments to __init__, this class can also generate instances of
    any subclass of SphericalParticle. Note that the additional args will not be Gaussian distributed, but fixed values.
    """
    def __init__(self,
                 number_of_particles: int,
                 position: Tuple[List[float], List[float]],
                 velocity: Tuple[List[float], List[float]],
                 radius: Tuple[float, float],
                 rho: float,
                 seed=None,
                 subclass: Type[SphericalParticle] = SphericalParticle,
                 **subclass_kwargs):
        self.particles = []
        if number_of_particles < 1:
            raise ValueError('The number of particles must be >= 1')

        self.number_of_particles = number_of_particles
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.radius = np.array(radius)
        self.rho = rho
        self.seed = seed
        self.subclass = subclass
        self.subclass_kwargs = subclass_kwargs

        self.number_of_dimensions = len(self.position[0])

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"Incompatible number of dimensions: {number_of_dimensions}, "
                f"the position mu/sigma description passed at construction was {self.number_of_dimensions}-dimensional."
            )

    def generate_particles(self, start_time: float = 0.0):
        """
        Given the parameters of the normal distribution this function generates
        the initial positions and momentum of the particles
        """
        if self.seed is not None:
            np.random.seed(self.seed)

        position = np.random.normal(self.position[0], self.position[1],
                                    (self.number_of_particles, self.number_of_dimensions))
        velocity = np.random.normal(self.velocity[0], self.velocity[1],
                                    (self.number_of_particles, self.number_of_dimensions))
        r = np.random.normal(self.radius[0], self.radius[1], self.number_of_particles)

        self.particles = []
        for i in range(self.number_of_particles):
            inst = self.subclass(
                identifier=i,
                position=np.concatenate([position[i], velocity[i]]),
                start_time=start_time,
                radius=r[i],
                rho=self.rho,
                **self.subclass_kwargs
            )
            self.particles.append(inst)

        return self.particles


class LinearlyDistributedSphericalSource(Source):
    def __init__(self,
                 number_of_particles: int,
                 position_intervals: List[Tuple[float, float]],
                 velocity_intervals: List[Tuple[float, float]],
                 radius_interval: Tuple[float, float],
                 rho: float,
                 seed=None,
                 radially_symmetrical: bool = False,
                 subclass: Type[SphericalParticle] = SphericalParticle,
                 **subclass_kwargs):
        self.particles = []
        if number_of_particles < 1:
            raise ValueError('The number of particles must be >= 1')

        self.number_of_particles = number_of_particles
        self.position_intervals = position_intervals
        self.velocity_intervals = velocity_intervals
        self.radius_interval = radius_interval
        self.rho = rho
        self.seed = seed
        self.radially_symmetrical = radially_symmetrical
        self.subclass = subclass
        self.subclass_kwargs = subclass_kwargs

        self.number_of_dimensions = len(self.position_intervals)

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"Incompatible number of dimensions: {number_of_dimensions}, "
                f"the position description passed at construction was {self.number_of_dimensions}-dimensional."
            )

    def _rotate(self, positions):
        size = positions.shape[1]
        rotated = (positions[0] + positions[1]*1j) * np.exp(np.random.random(size) * 2 * np.pi * 1j)
        return np.array([rotated.real, rotated.imag, positions[2]])

    def generate_particles(self, start_time: float = 0.0):
        """
        Given the parameters of the normal distribution this function generates
        the initial positions and momentum of the particles
        """
        if self.seed is not None:
            np.random.seed(self.seed)

        positions = np.array([
            np.linspace(pi[0], pi[1], self.number_of_particles) for pi in self.position_intervals
        ])
        if self.number_of_dimensions == 3 and self.radially_symmetrical:
            # Randomly rotate the positions around the Z axis for 3D radially symmetrical simulations
            positions = self._rotate(positions)
        positions = positions.transpose()

        velocities = np.array([
            np.linspace(vi[0], vi[1], self.number_of_particles) for vi in self.velocity_intervals
        ]).transpose()
        r = np.linspace(self.radius_interval[0], self.radius_interval[1], self.number_of_particles)

        self.particles = []
        for i in range(self.number_of_particles):
            inst = self.subclass(
                identifier=i,
                position=np.concatenate([positions[i], velocities[i]]),
                start_time=start_time,
                radius=r[i],
                rho=self.rho,
                **self.subclass_kwargs
            )
            self.particles.append(inst)

        return self.particles


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""