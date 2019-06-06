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

from typing import Tuple, Optional, List, Type, Dict, Any

import numpy as np

from cminject.experiment_refactored.definitions.base_classes import Particle, Boundary, Detector, Source, Calculator

infinite_interval = (float('-inf'), float('inf'))

# Can be used to describe the boundary of an object along an axis that has no extent along this axis.
empty_interval = (float('inf'), float('-inf'))


class TrajectoryCalculator(Calculator):
    def calculate(self, particle: Particle, time: float) -> None:
        particle.trajectory.append(
            np.concatenate([[time], particle.position])
        )


class DensityMapCalculator(Calculator):
    def __init__(self,
                 results: Dict[str, Any],
                 extent: np.array,
                 delta: float):
        bins = np.ceil(extent / delta).astype('int64')
        shape = bins[:, 1] - bins[:, 0]
        results['density_map'] = np.zeros(shape, dtype='int64')

        self.results = results
        self.extent = extent
        self.shape = shape
        self.delta = delta

        from scipy.interpolate import interp1d
        interpolators = []
        for i, axis in enumerate(extent):
            x = axis
            y = np.array([0, shape[i]-1])
            interp = interp1d(x, y, kind='linear')
            interpolators.append(interp)
        self.interpolators = interpolators

    def calculate(self, particle: Particle, time: float) -> None:
        bin_ = np.round([
            self.interpolators[i](particle.position[i])
            for i in range(len(particle.position))
        ]).astype(int)
        try:
            self.results['density_map'][bin_[0], bin_[1], bin_[2]] += 1
        except IndexError:
            pass


class SimpleZBoundary(Boundary):
    def __init__(self, z_minmax: Tuple[float, float]):
        z_min, z_max = z_minmax
        if z_min > z_max:
            raise ValueError("z_min must be < z_max!")
        self.z_min = z_min
        self.z_max = z_max
        self._z_boundary = (z_min, z_max)
        super().__init__()

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, particle: Particle):
        return self.z_min <= particle.position[2] <= self.z_max


class CuboidBoundary(SimpleZBoundary):
    def __init__(self, x_minmax: Tuple[float, float], y_minmax: Tuple[float, float], z_minmax: Tuple[float, float]):
        self.x_min, self.x_max = x_minmax
        self.y_min, self.y_max = y_minmax
        super().__init__(z_minmax)

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.x_min <= particle.position[0] <= self.x_max and\
               self.y_min <= particle.position[1] <= self.y_max and\
               self.z_min <= particle.position[2] <= self.z_max


class InfiniteBoundary(Boundary):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def is_particle_inside(self, particle: Particle) -> bool:
        return True


class SimpleZDetector(Detector):
    def __init__(self, identifier: int, z_position: float,):
        self.z_position = z_position
        self._z_boundary = (z_position, z_position)
        self.particle_distances = {}
        super().__init__(identifier=identifier)

    def _has_particle_reached_detector(self, particle: Particle) -> bool:
        reached = False
        prev_distance = self.particle_distances.get(particle.identifier, None)
        curr_distance = particle.position[2] - self.z_position
        if prev_distance and np.sign(curr_distance) != np.sign(prev_distance):
            reached = True
        self.particle_distances[particle.identifier] = curr_distance
        return reached

    def _hit_position(self, particle: Particle) -> Optional[np.array]:
        x, y, z, u, v, w = particle.position
        d = abs(z) - abs(self.z_position)
        t = d / abs(w)
        x = x - u * t
        y = y - v * t
        return np.array([x, y, self.z_position])

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary


class SphericalParticle(Particle):
    """
    A simple spherical particle that has a radius and a density (rho).
    """

    def __init__(self, *args, radius: float, rho: float, **kwargs):
        self.radius = radius
        self.rho = rho
        # self.mass is stored by the super() call.
        super().__init__(*args, **kwargs)

    def calculate_mass(self) -> float:
        return self.rho * 4 / 3 * np.pi * (self.radius ** 3)

    @property
    def properties(self) -> np.array:
        return np.concatenate([
            [
                self.identifier,
                self.time_of_flight,
                self.radius,
                self.mass
            ],
            self.initial_position,
        ])

    @property
    def properties_description(self) -> List[str]:
        return [
            'ID', 't', 'r', 'm',
            'x_i', 'y_i', 'z_i', 'u_i', 'v_i', 'w_i',
        ]


class ThermallyConductiveSphericalParticle(SphericalParticle):
    """
    Lacking a better name (so far), this class extends on the SphericalParticle with a thermal conductivity
    that is required by e.g. the photophoretic force and for thermal calculations.
    """
    def __init__(self, *args,
                 thermal_conductivity: float = 6.3,
                 temperature: float = 298.0,
                 **kwargs):
        self.thermal_conductivity = thermal_conductivity
        self.temperature = temperature
        self.collision_temperature = temperature

        self.time_to_liquid_n = float('inf')
        self.collision_time_to_liquid_n = float('inf')
        super().__init__(*args, **kwargs)

    @property
    def properties(self) -> np.array:
        basic_phase = super().properties
        return np.concatenate([
            basic_phase,
            [
                self.temperature,
                self.collision_temperature,
                self.time_to_liquid_n,
                self.collision_time_to_liquid_n,
            ]
        ])

    @property
    def properties_description(self) -> List[str]:
        return super().properties_description + [
            'T', 'T_c', 't_Ln', 't_Ln_c'
        ]


class GaussianSphericalSource(Source):
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
        super().__init__(position=position)

    def generate_particles(self, start_time: float = 0.0):
        """
        Given the parameters of the normal distribution this function generates
        the initial positions and momentum of the particles
        """
        if self.seed is not None:
            np.random.seed(self.seed)

        x, y, z = [
            np.random.normal(self.position[0][i], self.position[1][i], self.number_of_particles)
            for i in range(3)
        ]
        vx, vy, vz = [
            np.random.normal(self.velocity[0][i], self.velocity[1][i], self.number_of_particles)
            for i in range(3)
        ]
        r = np.random.normal(self.radius[0], self.radius[1], self.number_of_particles)

        self.particles = [
            self.subclass(
                identifier=i,
                position=np.array([x[i], y[i], z[i], vx[i], vy[i], vz[i]]),
                start_time=start_time,
                radius=r[i],
                rho=self.rho,
                **self.subclass_kwargs
            )
            for i in range(self.number_of_particles)
        ]
        return self.particles


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
