import numpy as np
from typing import Tuple
from cminject.definitions.base import Device
from cminject.definitions.fields.stark_field import B_StarkField
from cminject.definitions.boundaries import StarkBoundary
from cminject.definitions.particles import Molecule

class StarkDeflector(Device):

    def __init__(self, filename:str, z_minmax: Tuple[float, float], field_limit: float, *args, **kwargs):

        field = B_StarkField(filename=filename)
        boundary = StarkBoundary(field=field,z_minmax=z_minmax, field_limit=field_limit)
        super().__init__(fields=[field], boundary=boundary)
        # maybe I also need to create the rods of the field because they are also boundaries
    def calculate_acceleration(self, particle: Molecule, time: float):

        return self.fields[0].calculate_acceleration(particle, time)
