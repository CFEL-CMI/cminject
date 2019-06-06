from typing import Optional, Tuple

import numpy as np
from cminject.experiment_refactored.definitions.base import Detector, Particle


class SimpleZDetector(Detector):
    def __init__(self, identifier: int, z_position: float,):
        self.z_position = z_position
        self._z_boundary = (z_position, z_position)
        self.particle_distances = {}
        super().__init__(identifier=identifier)

    def _has_particle_reached_detector(self, particle: Particle) -> bool:
        reached = False
        prev_distance = self.particle_distances.get(particle.identifier, None)
        curr_distance = particle.position[2] - self.z_position
        if prev_distance and np.sign(curr_distance) != np.sign(prev_distance):
            reached = True
        self.particle_distances[particle.identifier] = curr_distance
        return reached

    def _hit_position(self, particle: Particle) -> Optional[np.array]:
        x, y, z, u, v, w = particle.position
        d = abs(z) - abs(self.z_position)
        t = d / abs(w)
        x = x - u * t
        y = y - v * t
        return np.array([x, y, self.z_position])

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary