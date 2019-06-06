from typing import Tuple, List, Type

import numpy as np
from cminject.experiment_refactored.definitions.base import Source
from cminject.experiment_refactored.definitions.particles import SphericalParticle


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