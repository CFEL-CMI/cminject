from typing import Tuple, List
import numpy as np

from cminject.experiment_refactored.definitions.base import Boundary, Particle, infinite_interval


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

    def is_particle_inside(self, particle: Particle):
        return self.z_min <= particle.position[self.number_of_dimensions] <= self.z_max


class CuboidBoundary(SimpleZBoundary):
    def __init__(self, intervals: List[Tuple[float, float]]):
        self.intervals = np.array(intervals)
        super().__init__(intervals[-1])  # last dimension is assumed to be Z

    def set_number_of_dimensions(self, number_of_dimensions: int):
        interval_dimensions = len(self.intervals)
        if number_of_dimensions != interval_dimensions:
            raise ValueError(
                f"Incompatible number of dimensions: {number_of_dimensions}, "
                f"the cuboid is defined within {interval_dimensions}-dimensional space"
            )
        self.number_of_dimensions = number_of_dimensions

    def is_particle_inside(self, particle: Particle) -> bool:
        pos = particle.position[:self.number_of_dimensions]
        return np.all(
            (self.intervals[:, 0] <= pos) * (pos <= self.intervals[:, 1])
        )


class InfiniteBoundary(Boundary):
    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def is_particle_inside(self, particle: Particle) -> bool:
        return True
