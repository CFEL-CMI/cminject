from typing import Tuple

from cminject.experiment_refactored.definitions.base import Boundary, Particle, infinite_interval


class SimpleZBoundary(Boundary):
    def __init__(self, z_minmax: Tuple[float, float]):
        z_min, z_max = z_minmax
        if z_min > z_max:
            raise ValueError("z_min must be < z_max!")
        self.z_min = z_min
        self.z_max = z_max
        self._z_boundary = (z_min, z_max)
        super().__init__()

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, particle: Particle):
        return self.z_min <= particle.position[2] <= self.z_max


class CuboidBoundary(SimpleZBoundary):
    def __init__(self, x_minmax: Tuple[float, float], y_minmax: Tuple[float, float], z_minmax: Tuple[float, float]):
        self.x_min, self.x_max = x_minmax
        self.y_min, self.y_max = y_minmax
        super().__init__(z_minmax)

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.x_min <= particle.position[0] <= self.x_max and\
               self.y_min <= particle.position[1] <= self.y_max and\
               self.z_min <= particle.position[2] <= self.z_max


class InfiniteBoundary(Boundary):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def is_particle_inside(self, particle: Particle) -> bool:
        return True