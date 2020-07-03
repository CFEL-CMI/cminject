#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of CMInject
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this program for scientific work, you should correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

from abc import ABC
from typing import List, Tuple, Any

import numpy as np

from cminject.definitions.base import ZBounded
from cminject.definitions.boundaries.base import Boundary
from cminject.definitions.fields.base import Field
from cminject.definitions.particles.base import Particle
from cminject.global_config import GlobalConfig, ConfigKey, ConfigSubscriber


class Device(ZBounded, ConfigSubscriber, ABC):
    """
    A combination of Field instances and a Boundary instance.
    Used to model real-world devices in an experiment setup.
    """
    def __init__(self, fields: List[Field], boundary: Boundary):
        """
        Constructor for Device.

        :param fields: The Fields that are present in this device.
        :param boundary: The Boundary this device has.
        """
        self.fields = fields
        self.boundary = boundary
        self.number_of_dimensions = None
        super().__init__()

        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            self.number_of_dimensions = value

    def calculate_acceleration(self, particle: Particle, time: float) -> np.array:
        # We've chosen to use this implementation as a default, as for n <= 10 fields, it is roughly 10 times faster
        # than using np.sum([f.calculate_acceleration(particle,time) for f in self.fields], axis=0). Devices with >10
        # fields are assumed to be rare.
        # This function should be overridden in a subclass if > 10 fields are used on this Device.
        # TODO reasons for this speed difference and numpy's worse performance should be investigated.
        acceleration = np.zeros(self.number_of_dimensions)
        for field in self.fields:
            acceleration += field.calculate_acceleration(particle, time)
        return acceleration

    def is_particle_inside(self, particle_position: np.array, time: float) -> bool:
        """
        Returns whether the passed Particle is inside this Device or not.
        Defaults to returning what the Boundary's method with the same name returns, but should be overridden
        if a more complex decision is required.

        :param particle_position: The particle's position.
        :param time: The current time.
        :return: True if the Particle inside this Device, False if it is not (or if this is unknown).
        """
        return self.boundary.is_particle_inside(particle_position, time)

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
