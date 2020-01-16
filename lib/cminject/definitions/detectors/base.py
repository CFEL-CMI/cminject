from abc import ABC, abstractmethod
from typing import Optional

import numpy as np

from cminject.definitions.base import NDimensional, ZBounded
from cminject.definitions.particles.base import Particle
from cminject.definitions.util import ParticleDetectorHit


class Detector(NDimensional, ZBounded, ABC):
    """
    Can tell whether a Particle has hit it, and _if_ it has hit it, can also tell where this occurred.

    As opposed to real detectors like fluorescent screens, these Detectors must never have an effect on particles.
    This allows us to, for instance, make an image of a particle distribution at different places of the experiment,
    without terminating the path of any particle.

    Detector subclasses can - and should - freely choose how to internally do the calculation of the hit position.
    """
    def __init__(self, identifier: int):
        """
        The constructor for Detector.

        :param identifier: A unique identifier for this detector. Used to find which hit occurred on which detector.
        """
        self.identifier = identifier
        super().__init__()

    def try_to_detect(self, particle: Particle) -> bool:
        """
        Tries to detect a particle.

        If the particle is considered having hit the detector ("could be detected")
        based on its current state, then a DetectorHit is stored on the particle's `detector_hits` property
        and True is returned. Otherwise, this function has no effect and returns False.

        :param particle: A Particle instance.
        :return: True if the particle has reached this detector, False otherwise.
        """
        if self._has_particle_reached_detector(particle.identifier, particle.position):
            hit_position = self._hit_position(particle.position)
            hit = ParticleDetectorHit(hit_position=hit_position, particle=particle)
            if self.identifier in particle.detector_hits:
                particle.detector_hits[self.identifier].append(hit)
            else:
                particle.detector_hits[self.identifier] = [hit]
            return True
        else:
            return False

    @abstractmethod
    def _has_particle_reached_detector(self, particle_identifier: int, position_velocity: np.array) -> bool:
        """Tells whether a Particle has reached this Detector. Must return True/False."""
        pass

    @abstractmethod
    def _hit_position(self, position_velocity: np.array) -> Optional[np.array]:
        """
        Can return a hit position on this detector for a Particle, but might
        also return None (if the Particle isn't considered having reached this detector).

        :param position_velocity: Position and velocity of the particle.
        :return: An (n,)-shaped numpy array describing the hit position if there is a hit position to calculate,
            None otherwise. n is the number of spatial dimensions in the experiment.
        """
        pass