from typing import Callable, Tuple

import numpy as np

from cminject.definitions import Field, Particle
from cminject.definitions.base import infinite_interval


class FunctionField(Field):
    """
    A simple stateless Field based on a function to calculate the acceleration.

    Implies infinite extent of the field, i.e. a particle is never considered to be outside the field's boundary.
    The number of dimensions of this field is not validated, meaning the given function must be able to handle
    particles of the experiment dimensionality.
    """
    def calculate_acceleration(self, particle: Particle, time: float) -> np.array:
        return self.function(particle, time)

    def __init__(self, function: Callable[[Particle, float], np.array]):
        self.function = function

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass
