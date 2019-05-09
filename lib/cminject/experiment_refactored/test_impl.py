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

from typing import List, Tuple, Optional

import numpy as np

from cminject.experiment_refactored.experiment import Experiment
from cminject.experiment_refactored.base_classes import\
    Particle, Source, Field, Boundary, Device, Detector

infinite_interval = (float('-inf'), float('inf'))


class SphericalParticle(Particle):
    def __init__(self, identifier, position, velocity, radius, rho):
        self.radius = radius
        self.rho = rho
        super().__init__(identifier, position, velocity)

    def calculate_mass(self) -> float:
        return self.rho * 4 / 3 * np.pi * (self.radius ** 3)


class OneToZeroBoundary(Boundary):
    def get_z_boundary(self) -> Tuple[float, float]:
        return 0.0, 1.0

    def is_particle_inside(self, particle: Particle):
        return 0.0 <= particle.position[2] <= 1.0


class GravityForceField(Field):
    def __init__(self):
        super().__init__()
        self.g = -9.81 * 20  # m/s^2

    def get_z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def calculate_acceleration(self,
                               position_and_velocity: np.array,
                               time: float,
                               particle: Particle = None) -> np.array:
        return np.array([
            0, 0, self.g
        ])


class ExampleDevice(Device):
    def __init__(self):
        super().__init__(field=GravityForceField(), boundary=OneToZeroBoundary())


class AtZEqualsOneHalfDetector(Detector):
    def has_reached_detector(self, particle: Particle) -> bool:
        if abs(particle.position[2] - 0.5) < 1e-3 and particle.position[2] >= 0.5:
            return True

    def get_hit_position(self, particle: Particle) -> Optional[np.array]:
        return particle.position

    def get_z_boundary(self) -> Tuple[float, float]:
        return 0.5, 0.5


class GaussianSphericalSource(Source):
    def __init__(self,
                 number_of_particles: int,
                 position: Tuple[List[float], List[float]],
                 velocity: Tuple[List[float], List[float]],
                 radius: Tuple[float, float] = (1.0e-6, 1.0e-7),
                 rho: float = 1.0,
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


def run_example_experiment(vz, do_profiling=False):
    devices = [ExampleDevice()]
    detectors = [AtZEqualsOneHalfDetector(identifier=1)]
    sources = [GaussianSphericalSource(
        10,
        position=([0.0, 0.0, 0.1], [0.0002, 0.0002, 0.0000]),
        velocity=([0.1, 0.0, vz], [0.0002, 0.0002, 0.00000]),
        radius=(0.3, 0.2),
        seed=100
    )]
    experiment = Experiment(
        devices=devices,
        detectors=detectors,
        sources=sources,
        t_end=0.1,
        dt=0.00001,
        track_trajectories=True
    )
    result = experiment.run(do_profiling=do_profiling)
    return result


def plot_trajectories(experiment_result):
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

    xs0, ys0, zs0 = [[p.initial_position[i] for p in result_list] for i in range(3)]
    xs, ys, zs = [[p.position[i] for p in result_list] for i in range(3)]
    r = [50 * p.radius for p in result_list]
    color = [
        'red' if p.lost else 'green'
        for p in result_list
    ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.scatter(xs0, ys0, zs=zs0, s=r)
    plt.scatter(xs, ys, zs=zs, s=r, c=color)

    for p in result_list:
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
            plt.scatter(hits[:, 0], hits[:, 1], zs=hits[:, 2], s=10, color='blue')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()


if __name__ == '__main__':
    import sys
    do_profiling = (len(sys.argv) > 1 and sys.argv[1] == 'profile')
    result_list = run_example_experiment(vz=13.0, do_profiling=do_profiling)
    plot_trajectories(result_list)


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""