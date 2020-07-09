#!/usr/bin/env python3
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

from abc import ABC, abstractmethod
from typing import Optional

import numpy as np

from cminject.definitions.base import ZBounded
from cminject.definitions.particles.base import Particle
from cminject.definitions.util import ParticleDetectorHit


class Detector(ZBounded, ABC):
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
            hit_position = self._hit_position(particle.phase_space_position)
            hit = ParticleDetectorHit(hit_position=hit_position, particle=particle)
            if self.identifier in particle.detector_hits:
                particle.detector_hits[self.identifier].append(hit)
            else:
                particle.detector_hits[self.identifier] = [hit]
            return True
        else:
            return False

    @abstractmethod
    def _has_particle_reached_detector(self, particle_identifier: int, position: np.array) -> bool:
        """Tells whether a Particle has reached this Detector. Must return True/False."""
        pass

    @abstractmethod
    def _hit_position(self, phase_space_position: np.array) -> Optional[np.array]:
        """
        Can return a hit position on this detector for a Particle, but might
        also return None (if the Particle isn't considered having reached this detector).

        :param phase_space_position: Position and velocity of the particle.
        :return: An (n,)-shaped numpy array describing the (spatial) hit position if there is a hit position
            to calculate, None otherwise. n is the number of spatial dimensions in the experiment.
        """
        pass

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
