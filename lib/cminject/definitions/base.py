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

__all__ = ['NDimensional', 'ZBounded']

from abc import ABC, abstractmethod
from typing import Tuple


class NDimensional(ABC):
    """
    A mixin for objects that are N-dimensional, or at least (may) act different depending on the number of spatial
    dimensions of the simulation setup they are used within.

    This method is abstract to force every relevant object in the experimental setup to store and respond appropriately
    to the dimensionality of the whole simulation. If your concrete implementation for any class that's derived from
    this does not care about dimensionality at all, simply implement `set_number_of_dimensions` and have it do nothing,
    i.e. just `pass` as an implementation.
    """
    @abstractmethod
    def set_number_of_dimensions(self, number_of_dimensions: int):
        """
        A method to change this object's data/state/implementation/... to conform to the choice of the simulation
        setup's number of spatial dimensions.

        What this means concretely really depends on the object type and implementation, but for example,
        after a `set_number_of_dimensions(2)` call:

        - A Source will now generate particles with a 2D distribution instead of a 3D distribution
        - A Field will now calculate a 2D acceleration instead of a 3D acceleration

        You can expect this method to be only called once, when the experiment is preparing to be run.
        From then on, the dimensionality of the object should be considered fixed.

        .. note::
            If your object can only handle certain numbers of dimensions (e.g. only 3D, or only 2D and 3D), it is
            good practice to raise a ValueError exception when this method is called with a dimensionality value
            your object can not work with. This will make it clear why the simulation setup does not make sense as
            it is, and how it must be changed before it can be run.

        :param number_of_dimensions: The number of spatial dimensions of the setup this object should work within.
        :return: Nothing.
        """
        pass


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
