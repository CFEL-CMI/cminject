from typing import Tuple

import numpy as np

from .base import Boundary
from cminject.definitions.fields.regular_grid_interpolation import RegularGridInterpolationField


class GridFieldBasedBoundary(Boundary):
    def set_number_of_dimensions(self, number_of_dimensions: int):
        d = number_of_dimensions
        fd = self.field.number_of_dimensions
        if d != fd:
            raise ValueError(f"Dimensionality {d} does not match associated field's dimensionality {fd}!")

    def __init__(self, field: RegularGridInterpolationField):
        self.field = field

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.field.z_boundary

    def is_particle_inside(self, position: np.array, time: float) -> bool:
        return self.field.is_particle_inside(position, time)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End: