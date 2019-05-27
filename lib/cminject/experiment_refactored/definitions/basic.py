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

from typing import Tuple, Optional, List, Type

import numpy as np

from cminject.experiment_refactored.definitions.base_classes import Particle, Boundary, Detector, Source


infinite_interval = (float('-inf'), float('inf'))


class SimpleZBoundary(Boundary):
    def __init__(self, z_min, z_max):
        if z_min > z_max:
            raise ValueError("z_min must be < z_max!")
        self.z_min = z_min
        self.z_max = z_max
        super().__init__()

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.z_min, self.z_max

    def is_particle_inside(self, particle: Particle):
        return self.z_min <= particle.position[2] <= self.z_max


class CuboidBoundary(Boundary):
    def __init__(self, x_minmax: Tuple[float, float], y_minmax: Tuple[float, float], z_minmax: Tuple[float, float]):
        self.x_min, self.x_max = x_minmax
        self.y_min, self.y_max = y_minmax
        self.z_min, self.z_max = z_minmax
        super().__init__()

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.z_min, self.z_max

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.x_min <= particle.position[0] <= self.x_max and\
               self.y_min <= particle.position[1] <= self.y_max and\
               self.z_min <= particle.position[2] <= self.z_max


class InfiniteBoundary(Boundary):
    def get_z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def is_particle_inside(self, particle: Particle) -> bool:
        return True


class SimpleZDetector(Detector):
    def __init__(self, identifier: int, z_position: float,):
        super().__init__(identifier=identifier)
        self.z_position = z_position
        self.particle_distances = {}

    def has_reached_detector(self, particle: Particle) -> bool:
        reached = False
        prev_distance = self.particle_distances.get(particle.identifier, None)
        curr_distance = particle.position[2] - self.z_position
        if prev_distance and np.sign(curr_distance) != np.sign(prev_distance):
            reached = True
        self.particle_distances[particle.identifier] = curr_distance
        return reached

    def get_hit_position(self, particle: Particle) -> Optional[np.array]:
        return particle.position

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.z_position, self.z_position


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


class ThermallyConductiveSphericalParticle(SphericalParticle):
    """
    Lacking a better name (so far), this class extends on the SphericalParticle with a thermal conductivity
    that is required by e.g. the photophoretic force.
    """
    def __init__(self, *args,
                 thermal_conductivity: float = 6.3,
                 temperature: float = 298.0,
                 **kwargs):
        self.thermal_conductivity = thermal_conductivity
        self.temperature = temperature
        self.collision_temperature = temperature

        self.time_to_liquid_n = None
        self.collision_time_to_liquid_n = None
        super().__init__(*args, **kwargs)


class GaussianSphericalSource(Source):
    """
    A Source generating a Gaussian distribution of SphericalParticles wrt position, velocity and radius.

    By passing `subclass` and any additional arguments to __init__, this class can also generate instances of
    any subclass of SphericalParticle. The additional args will not be Gaussian distributed, but fixed values.
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
                position=np.array([x[i], y[i], z[i]]),
                velocity=np.array([vx[i], vy[i], vz[i]]),
                start_time=start_time,
                radius=r[i],
                rho=self.rho,
                **self.subclass_kwargs
            )
            for i in range(self.number_of_particles)
        ]
        return self.particles


def plot_particles(experiment_result, plot_trajectories=False):
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    xs0, ys0, zs0 = [[p.initial_position[i] for p in experiment_result] for i in range(3)]
    xs, ys, zs = [[p.position[i] for p in experiment_result] for i in range(3)]
    color = [
        'red' if p.lost and not p.reached else 'green'
        for p in experiment_result
    ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.scatter(xs, ys, zs=zs, s=3, c=color)
    plt.scatter(xs0, ys0, zs=zs0, s=3, c='black')

    for p in experiment_result:
        for detector_id, hits in p.detector_hits.items():
            hits = np.array(hits)
            plt.scatter(hits[:, 0], hits[:, 1], zs=hits[:, 2], s=10, color='yellow')

        if plot_trajectories and p.trajectory:
            trajectory = np.array(p.trajectory)
            ts = trajectory[:, 0]
            ps = trajectory[:, 1:4]
            xs = ps[:, 0]
            ys = ps[:, 1]
            zs = ps[:, 2]

            vs = trajectory[:, 4:]
            plt.plot(xs, ys, zs=zs, color='grey')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
