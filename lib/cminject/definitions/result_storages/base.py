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
from typing import List, Iterable, Dict, Optional

from cminject.definitions.particles.base import Particle


class ResultStorage(ABC):
    """
    An object to store the results of an experiment in some fashion. MUST implement store_results, and MAY implement
    convenience methods to read from the storage again (e.g. a method to get all particle trajectories from a file).
    """
    @abstractmethod
    def store_results(self, particles: List[Particle]):
        """
        Stores the results of an experiment (which are always a list of modified Particle instances).

        :param particles: The list of particles, each in the state of after running a simulation.
        """
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    @abstractmethod
    def get_dimensions(self) -> int:
        pass

    @abstractmethod
    def get_identifiers(self) -> Iterable[str]:
        pass

    @abstractmethod
    def get_initial_positions(self) -> Optional[Iterable]:
        pass

    @abstractmethod
    def get_final_positions(self) -> Optional[Iterable]:
        pass

    @abstractmethod
    def get_trajectories(self) -> Optional[Iterable]:
        pass

    @abstractmethod
    def get_detectors(self) -> Optional[Dict[str, Iterable]]:
        pass

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
