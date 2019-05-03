from typing import List, Tuple, Optional

import numpy as np

from cminject.experiment_refactored.experiment import Experiment
from cminject.experiment_refactored.base_classes import\
    Particle, Source, Field, Boundary, Device, Detector

infinite_interval = (float('-inf'), float('inf'))


class SphericalParticle(Particle):
    def __init__(self, identifier, position, velocity, radius, rho):
        mass = rho * 4/3 * np.pi * (radius**3)
        super().__init__(identifier, position, velocity, mass)

        self.radius = radius
        self.rho = rho


class OneToZeroBoundary(Boundary):
    def is_particle_inside(self, particle: Particle):
        z = particle.position[2]
        return 0.0 <= z <= 1.0

    def get_z_boundary(self) -> Tuple[float, float]:
        return 0.0, 1.0


class GravityField(Field):
    def __init__(self):
        super().__init__(needs_full_particle_data=True)
        self.g = -9.81  # m/s^2

    def get_z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def calculate_acceleration(self,
                               position_and_velocity: np.array,
                               particle: Particle = None) -> np.array:
        return np.array([0, 0, self.g * particle.mass])  # TODO this is physical nonsense, just for demo


class GravityDevice(Device):
    def __init__(self):
        super().__init__(field=GravityField(), boundary=OneToZeroBoundary())


class BlindDetector(Detector):
    def has_reached_detector(self, particle: Particle) -> bool:
        return False

    def get_hit_position(self, particle: Particle) -> Optional[np.array]:
        return None

    def get_z_boundary(self) -> Tuple[float, float]:
        return infinite_interval


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


def run_gravity_experiment(vz):
    devices = [GravityDevice()]
    detectors = [BlindDetector()]
    sources = [GaussianSphericalSource(
        100,
        position=([0.0, 0.0, 0.0], [0.0002, 0.0002, 0.0000]),
        velocity=([0.0, 0.0, vz], [0.001, 0.001, 0.00002]),
        radius=(0.1, 0.2)
    )]
    experiment = Experiment(
        devices=devices,
        detectors=detectors,
        sources=sources,
        t_end=0.1
    )
    return experiment.run()


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

    result_list = run_gravity_experiment(0.0)
    xs0, ys0, zs0 = [[p.initial_position[i] for p in result_list] for i in range(3)]
    xs, ys, zs = [[p.position[i] for p in result_list] for i in range(3)]
    r = [p.mass * 5 for p in result_list]

    if True:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.scatter(xs0, ys0, zs=zs0, s=r)
        plt.scatter(xs, ys, zs=zs, s=r)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()
