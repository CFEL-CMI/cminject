from abc import ABC, abstractmethod
from typing import List

from cminject.definitions.particles.base import Particle


class ResultStorage(ABC):
    """
    An object to store the results of an experiment in some fashion. MUST implement store_results, and MAY implement
    convenience methods to read from the storage again (e.g. a method to get all particle trajectories from a file).
    """
    @abstractmethod
    def store_results(self, particles: List[Particle]):
        """
        Stores the results of an experiment (which are always a list of modified Particle instances).

        :param particles: The list of particles, each in the state of after running a simulation.
        """
        pass