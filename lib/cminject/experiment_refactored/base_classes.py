from abc import ABC, abstractmethod
from typing import Tuple, List, Any, Optional, Dict, Set

import numpy as np


class Particle(ABC):
    def __init__(self, identifier: Any,
                 position: np.array, velocity: np.array,
                 mass: float):
        self.identifier = identifier
        self.lost = False
        self.detector_hits: Dict[int, List[np.array]] = {}
        self.position = position
        self.initial_position = np.copy(position)
        self.velocity = velocity
        self.mass = mass
        self.trajectory = []

    def store_hit(self, identifier: int, position: np.array):
        if identifier in self.detector_hits:
            self.detector_hits[identifier].append(position)
        else:
            self.detector_hits[identifier] = [position]
        pass


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
    def __init__(self, identifier: int):
        self.identifier = identifier

    @abstractmethod
    def has_reached_detector(self, particle: Particle) -> bool:
        pass

    @abstractmethod
    def get_hit_position(self, particle: Particle) -> Optional[np.array]:
        pass


class Boundary(ZBoundedMixin, ABC):
    @abstractmethod
    def is_particle_inside(self, particle: Particle) -> bool:
        pass


class Field(ZBoundedMixin, ABC):
    @abstractmethod
    def calculate_acceleration(self,
                               position_and_velocity: np.array,
                               particle: Particle,
                               time: float) -> np.array:
        pass


class Device(ZBoundedMixin, ABC):
    def __init__(self, field: Field, boundary: Boundary):
        self.field = field
        self.boundary = boundary

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.boundary.get_z_boundary()

    def is_particle_inside(self, particle: Particle):
        return self.boundary.is_particle_inside(particle)
