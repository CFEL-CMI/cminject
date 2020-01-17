from abc import ABC, abstractmethod

import numpy as np

from cminject.definitions.base import NDimensional, ZBounded
from cminject.definitions.particles.base import Particle


class Field(NDimensional, ZBounded, ABC):
    """
    A Field interacting with Particle objects.
    """
    @abstractmethod
    def calculate_acceleration(self,
                               particle: Particle,
                               time: float) -> np.array:
        """
        Calculates an acceleration for one particle based on the particle's current properties and the current time.
        This acceleration will be integrated for in each time step and thus "applied" to the particle.

        :param particle: The Particle to calculate this Field's acceleration for.
        :param time: The time to calculate the acceleration for.
        :return: A (n,)-shaped numpy array describing the acceleration exerted on the particle. n is the number of
            spatial dimensions of the experiment.
        """
        pass