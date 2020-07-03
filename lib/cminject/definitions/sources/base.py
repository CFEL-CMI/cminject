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
from typing import List

from cminject.definitions.particles.base import Particle


class Source(ABC):
    """
    A source of particles, to generate an initial position (and property) distribution.

    One can implement sources following different distributions and generating different particle types by subclassing
    this class and implementing generate_particles however they want.

    .. warning::
      Any implementation of Source should only generate exactly one type of particle, i.e. its generate_particles
      method should never return a list of particles where the (most specific) type of any two particles is different.
    """
    @abstractmethod
    def generate_particles(self, start_time: float = 0.0) -> List[Particle]:
        """
        Generates a list of particles. How this is done is entirely up to the subclass, by (this is mandatory!)
        implementing this method in some way.

        :return: A list of Particle instances.
        """
        pass

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
