from typing import Optional, Tuple, List

import numpy as np
from cminject.experiment_refactored.definitions.base import Detector, Particle


class SimpleZDetector(Detector):
    def __init__(self, identifier: int, z_position: float):
        self.z_position = z_position
        self._z_boundary = (z_position, z_position)
        self.particle_distances = {}

        self.number_of_dimensions = 0
        self._position_descriptions = []

        super().__init__(identifier=identifier)

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions == 3:
            self._hit_position = self._hit_position3
            self._position_descriptions = ['x_d', 'y_d', 'z_d']
        elif number_of_dimensions == 2:
            self._hit_position = self._hit_position2
            self._position_descriptions = ['r_d', 'z_d']
        elif number_of_dimensions == 1:
            self._hit_position = self._hit_position1
            self._position_descriptions = ['z_d']
        else:
            raise ValueError("Dimensionality has to be in [1, 2, 3].")

        self.number_of_dimensions = number_of_dimensions

    def position_descriptions(self) -> List[str]:
        return self._position_descriptions

    def _has_particle_reached_detector(self, particle: Particle) -> bool:
        reached = False
        prev_distance = self.particle_distances.get(particle.identifier, None)
        curr_distance = particle.position[self.number_of_dimensions - 1] - self.z_position
        if prev_distance and np.sign(curr_distance) != np.sign(prev_distance):
            reached = True
        self.particle_distances[particle.identifier] = curr_distance
        return reached

    def _hit_position(self, particle: Particle) -> Optional[np.array]:
        raise Exception("This should never be called, the implementation should be swapped out by the constructor!")

    def _hit_position3(self, particle: Particle) -> Optional[np.array]:
        x, y, z, u, v, w = particle.position
        d = abs(z) - abs(self.z_position)
        t = d / abs(w)
        x = x - u * t
        y = y - v * t
        return np.array([x, y, self.z_position])

    def _hit_position2(self, particle: Particle) -> Optional[np.array]:
        r, z, vr, vz = particle.position
        d = abs(z) - abs(self.z_position)
        t = d / abs(vz)
        r = r - vr * t
        return np.array([r, self.z_position])

    def _hit_position1(self, particle: Particle) -> Optional[np.array]:
        return np.array([self.z_position])

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary