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

from typing import Tuple, Optional, List

import numpy as np

from .base_classes import Particle, Boundary, Detector, Source


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


class SimpleZDetector(Detector):
    def __init__(self, identifier: int, z_position: float, epsilon: float = 1e-3):
        super().__init__(identifier=identifier)
        self.z_position = z_position
        self.epsilon = epsilon

    def has_reached_detector(self, particle: Particle) -> bool:
        if abs(particle.position[2] - self.z_position) < self.epsilon\
                and particle.position[2] >= self.z_position:
            return True

    def get_hit_position(self, particle: Particle) -> Optional[np.array]:
        return particle.position

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.z_position, self.z_position


class SphericalParticle(Particle):
    def __init__(self, identifier, position, velocity, radius, rho):
        self.radius = radius
        self.rho = rho
        super().__init__(identifier, position, velocity)

    def calculate_mass(self) -> float:
        return self.rho * 4 / 3 * np.pi * (self.radius ** 3)


class GaussianSphericalSource(Source):
    def __init__(self,
                 number_of_particles: int,
                 position: Tuple[List[float], List[float]],
                 velocity: Tuple[List[float], List[float]],
                 radius: Tuple[float, float],
                 rho: float,
                 seed=None):
        self.particles = []
        if number_of_particles < 1:
            raise ValueError('The number of particles must be >= 1')

        self.number_of_particles = number_of_particles
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.radius = np.array(radius)
        self.rho = rho
        self.seed = seed
        super().__init__(position=position)

    def generate_particles(self):
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
            SphericalParticle(
                i,
                position=np.array([x[i], y[i], z[i]]),
                velocity=np.array([vx[i], vy[i], vz[i]]),
                radius=r[i],
                rho=self.rho
            )
            for i in range(self.number_of_particles)
        ]
        return self.particles


def plot_trajectories(experiment_result):
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

    xs0, ys0, zs0 = [[p.initial_position[i] for p in experiment_result] for i in range(3)]
    xs, ys, zs = [[p.position[i] for p in experiment_result] for i in range(3)]
    r = [50 * p.radius for p in experiment_result]
    color = [
        'red' if p.lost else 'green'
        for p in experiment_result
    ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.scatter(xs0, ys0, zs=zs0, s=r)
    plt.scatter(xs, ys, zs=zs, s=r, c=color)

    for p in experiment_result:
        if not p.trajectory:
            continue

        traj = np.array(p.trajectory)
        ts = traj[:, 0]
        ps = traj[:, 1:4]
        xs = ps[:, 0]
        ys = ps[:, 1]
        zs = ps[:, 2]

        vs = traj[:, 4:]
        plt.plot(xs, ys, zs=zs, color='grey')

        for detector_id, hits in p.detector_hits.items():
            hits = np.array(hits)
            plt.scatter(hits[:, 0], hits[:, 1], zs=hits[:, 2], s=10, color=
                        ['yellow', 'orange', 'red'][(detector_id - 1) % 3])

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