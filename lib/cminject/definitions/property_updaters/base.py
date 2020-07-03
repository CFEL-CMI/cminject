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

from abc import ABC, abstractmethod

from cminject.definitions.particles.base import Particle


class PropertyUpdater(ABC):
    """
    An object to update a Particle's properties based on its current state. Can do arbitrary calculations to determine
    the values the properties should be set to.

    .. warning::
        If the properties that a .update() call changes include the .position attribute, the update method
        MUST return True, and it SHOULD return False otherwise. If your updater doesn't adhere to this, expect these
        updated positions to have no effect, since the integrator will overwrite them.

    Should be used for all calculations of additional quantities beyond calculating an acceleration, so,
    beyond what a Field returns. For example, TrajectoryPropertyUpdater is a simple PropertyUpdater that just takes the
    position and appends it to the `trajectory` property on the particle.

    The `update` method is guaranteed by `Experiment` to be called exactly once per time step (as long as the
    particle is still being simulated).

    If you need access to properties of fields, detectors, etc., override the __init__ method and store a reference
    to each relevant object at construction of the PropertyUpdater instance.

    .. note::
        Taking care to avoid data clashes (different property updaters writing on the same property on a particle,
        overwriting each others' results) is up to the user / implementer. This cannot be avoided in a general way.
        If you know that another PropertyUpdater instance is writing to a specific particle property, don't write to it
        too, unless you know the order of execution of the PropertyUpdaters and are doing this completely intentionally.
    """
    @abstractmethod
    def update(self, particle: Particle, time: float) -> bool:
        """
        Updates a property of some Particle instance in some way.

        :param particle: The Particle instance.
        :param time: The current time.
        :return: True if the position of a particle was changed
        """
        pass

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
