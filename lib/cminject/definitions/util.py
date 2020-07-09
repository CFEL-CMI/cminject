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
Smaller utility definitions that have no place anywhere else.
"""
import numpy as np
from cminject.definitions.particles.base import Particle

"""
An interval with no extent. Should be used to describe the boundary of an object
along an axis that has no extent along this axis.

The implementation is [∞, -∞], which is actually an impossible interval.
There is a technical reason to this: We gather the total Z boundary of the simulation by getting the minimum and maximum
extent of all objects in Z direction, min() is called on all left values and max() is called on all right values.
Since ∞ is the neutral element of min() and -∞ is the neutral element of max(), this works out just as intended:
The object without extent has no effect on the calculation of the total Z boundary. 
"""
empty_interval = (float('inf'), float('-inf'))

"An interval with infinite extent."
infinite_interval = (float('-inf'), float('inf'))


class ParticleDetectorHit(object):
    """
    A small data container class, to store required information about a detection event of a particle by a detector.

    Contains :attr:`hit_state`, which is a full description of the particle's state, with the 'position' part
    overwritten by the position the detector detected the particle at.
    """
    def __init__(self, hit_position: np.array, particle: Particle):
        """
        The constructor for a ParticleDetectorHit.

        :param hit_position: An (n,)-shaped np.array representing the hit position of the particle on the detector.
            n is the number of spatial dimensions and must match that of the particle.
        :param particle: The particle which was detected.
        """
        dim_pos = hit_position.size
        dim_par = particle.number_of_dimensions
        if dim_pos != dim_par:
            raise ValueError(
                f"Hit position dimensionality ({dim_pos}) does not match the number of "
                f"dimensions the particle is simulated in ({dim_par})!"
            )

        hit_state = particle.as_array('all')
        hit_state['position'] = hit_position
        self._hit_state = hit_state

    @property
    def hit_state(self) -> np.array:
        """
        The state of the particle in the moment where the detector detected the particle.
        The 'position' part of this state refers to the (approximate) position of detection as returned by the detector.

        :return: The position of the hit.
        """
        return self._hit_state

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
