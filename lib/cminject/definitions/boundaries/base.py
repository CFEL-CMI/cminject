#!/usr/bin/env python
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

from abc import ABC, abstractmethod

import numpy as np

from cminject.definitions.base import NDimensional, ZBounded


class Boundary(NDimensional, ZBounded, ABC):
    """
    Can tell whether a Particle is inside of it or not.

    Making the total set of Boundary objects in an experiment closed is - for now - the job of the person implementing
    the experiment setup.
    """
    @abstractmethod
    def is_particle_inside(self, position: np.array, time: float) -> bool:
        """
        Tells whether the passed Particle is inside of this Boundary or not.

        :param position: The Particle's position to tell this for.
        :param time: The time to tell this for.
        :return: True if the particle is definitely inside this Boundary,
            False if it is not inside this Boundary or if this cannot be known.
        """
        pass


