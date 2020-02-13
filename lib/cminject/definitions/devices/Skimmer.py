import numpy as np
from cminject.definitions.base import Device
from cminject.definitions.boundaries import Ring

class Skimmer(Device):

    def __init__(self, radius: float, position: np.array, *args, **kwargs):

        boundary = Ring(radius=radius,position=position)
        super().__init__(fields=[None], boundary=boundary)
        # maybe I also need to create the rods of the field because they are also boundaries
    def calculate_acceleration(self, particle: Molecule, time: float):

        return 0
