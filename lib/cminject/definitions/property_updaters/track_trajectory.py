import numpy as np

from cminject.definitions.particles.base import Particle
from .base import PropertyUpdater


class TrajectoryPropertyUpdater(PropertyUpdater):
    """
    A simple property updater to append to the trajectory of the particle after each time step.
    """
    def update(self, particle: Particle, time: float) -> bool:
        particle.trajectory.append(
            np.concatenate([[time], particle.position])
        )
        return False

    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass