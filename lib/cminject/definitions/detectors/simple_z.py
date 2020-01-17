from typing import Optional, Dict, Any, Tuple

import numpy as np

from .base import Detector


class SimpleZDetector(Detector):
    def _hit_position(self, position_velocity: np.array) -> Optional[np.array]:
        raise Exception("This should never be called, the implementation should be swapped out by the constructor!")

    def __init__(self, identifier: int, z_position: float):
        self.z_position = z_position
        self._z_boundary = (z_position, z_position)
        self.particle_distances: Dict[Any, float] = {}

        self.number_of_dimensions = 0

        super().__init__(identifier=identifier)

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions == 3:
            self._hit_position = self._hit_position3
        elif number_of_dimensions == 2:
            self._hit_position = self._hit_position2
        else:
            raise ValueError(f"{self.__class__.__name__}'s Dimensionality has to be in [2, 3].")

        self.number_of_dimensions = number_of_dimensions

    def _has_particle_reached_detector(self, particle_identifier: int, position_velocity: np.array) -> bool:
        if position_velocity[self.number_of_dimensions - 1] == self.z_position:
            return True

        reached = False
        prev_distance = self.particle_distances.get(particle_identifier, None)
        curr_distance = position_velocity[self.number_of_dimensions - 1] - self.z_position
        if prev_distance and np.sign(curr_distance) != np.sign(prev_distance):
            reached = True
        self.particle_distances[particle_identifier] = curr_distance
        return reached

    def _hit_position3(self, position_velocity: np.array) -> Optional[np.array]:
        if position_velocity[self.number_of_dimensions-1] == self.z_position:
            return position_velocity

        x, y, z, u, v, w = position_velocity
        d = abs(z) - abs(self.z_position)
        t = d / abs(w)
        x = x - u * t
        y = y - v * t
        return np.array([x, y, self.z_position, u, v, w])

    def _hit_position2(self, position_velocity: np.array) -> Optional[np.array]:
        if position_velocity[self.number_of_dimensions - 1] == self.z_position:
            return position_velocity

        r, z, vr, vz = position_velocity
        d = abs(z) - abs(self.z_position)
        t = d / abs(vz)
        r = r - vr * t
        return np.array([r, self.z_position, vr, vz])

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary