from typing import Tuple

import numpy as np

from .base import Boundary


class SimpleZBoundary(Boundary):
    def __init__(self, z_minmax: Tuple[float, float]):
        z_min, z_max = z_minmax
        if z_min > z_max:
            raise ValueError("z_min must be < z_max!")
        self.z_min = z_min
        self.z_max = z_max
        self._z_boundary = (z_min, z_max)

        self.number_of_dimensions = 3
        super().__init__()

    def set_number_of_dimensions(self, number_of_dimensions: int):
        self.number_of_dimensions = number_of_dimensions

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, position: np.array, time: float):
        return self.z_min <= position <= self.z_max