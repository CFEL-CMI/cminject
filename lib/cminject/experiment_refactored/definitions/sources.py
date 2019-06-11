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

        self.number_of_dimensions = len(self.position[0])
        super().__init__(position=position)

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"Incompatible number of dimensions: {number_of_dimensions},"
                f"the position mu/sigma description passed at construction is {self.number_of_dimensions}-dimensional."
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
