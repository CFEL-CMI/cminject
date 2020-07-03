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

"""
The most fundamental of definitions, for deriving other base classes from.
"""

__all__ = ['ZBounded']

from abc import ABC, abstractmethod
from typing import Tuple

class ZBounded(ABC):
    """
    A mixin for objects in an experiment setup that are bounded in the Z direction.
    Classes deriving from this mixin will have to implement the get_z_boundary method.

    Objects implementing this mixin should be used to determine the minimum/maximum extent of the entire experiment
    setup in the Z direction, and then allows an Experiment to make a fast preliminary calculation about whether it
    can consider a particle lost based on this.
    """
    @property
    @abstractmethod
    def z_boundary(self) -> Tuple[float, float]:
        """
        Returns the Z boundary of this Z-bounded object.

        :return: A tuple of floats, the first entry being z_min, the second being z_max.
        """
        pass

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
