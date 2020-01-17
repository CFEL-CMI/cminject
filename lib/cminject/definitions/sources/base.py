from abc import ABC, abstractmethod
from typing import List

from cminject.definitions.base import NDimensional
from cminject.definitions.particles.base import Particle


class Source(NDimensional, ABC):
    """
    A source of particles, to generate an initial position (and property) distribution.

    One can implement sources following different distributions and generating different particle types by subclassing
    this class and implementing generate_particles however they want.
    """
    @abstractmethod
    def generate_particles(self, start_time: float = 0.0) -> List[Particle]:
        """
        Generates a list of particles. How this is done is entirely up to the subclass, by (this is mandatory!)
        implementing this method in some way.

        :return: A list of Particle instances.
        """
        pass