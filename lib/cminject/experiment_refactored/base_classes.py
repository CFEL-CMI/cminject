from abc import ABC, abstractmethod
from typing import Tuple, List, Any, Optional

import numpy as np


class Particle(ABC):
    def __init__(self, identifier: Any,
                 position: np.array, velocity: np.array,
                 mass: float):
        self.identifier = identifier
        self.position = position
        self.initial_position = np.copy(position)
        self.velocity = velocity
        self.mass = mass


class Source(ABC):
    def __init__(self, position: np.array):
        self.position = position

    @abstractmethod
    def generate_particles(self) -> List[Particle]:
        pass


class ZBoundedMixin:
    @abstractmethod
    def get_z_boundary(self) -> Tuple[float, float]:
        pass


class Detector(ZBoundedMixin, ABC):
    def __init__(self, interpolate_hit_position=False):
        self.interpolate_hit_position = interpolate_hit_position

    @abstractmethod
    def has_reached_detector(self, particle: Particle) -> bool:
        pass

    @abstractmethod
    def get_hit_position(self, particle: Particle) -> Optional[np.array]:
        pass


class Boundary(ABC, ZBoundedMixin):
    @abstractmethod
    def is_particle_inside(self, particle: Particle):
        pass


class Field(ZBoundedMixin, ABC):
    def __init__(self, needs_full_particle_data=False):
        self.needs_full_particle_data = needs_full_particle_data

    @abstractmethod
    def calculate_acceleration(self,
                               position_and_velocity: np.array,
                               particle: Particle) -> np.array:
        pass


class Device(ZBoundedMixin, ABC):
    def __init__(self, field: Field, boundary: Boundary):
        self.field = field
        self.boundary = boundary

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.boundary.get_z_boundary()

    def is_particle_inside(self, particle: Particle):
        return self.boundary.is_particle_inside(particle)
