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
from typing import Any, List, Dict
from functools import cached_property

import numpy as np

from cminject.global_config import ConfigSubscriber, ConfigKey, GlobalConfig


class Particle(ConfigSubscriber, ABC):
    """
    Describes a particle whose trajectory we want to simulate.
    It is first and foremost a data container, and it and its subclasses should be written and used as such.

    It can be read by any part of the code, but should only be written to by instances of Experiment
    (or subclasses thereof).

    This class declares a few basic properties within __init__ that we expect every particle
    to have:
    - identifier (a unique id),
    - lost (a flag storing whether the particle is considered lost)
    - detector_hits (a dict mapping detector identifiers to hit lists)
    # TODO more props
    - position (the particle's position)
    - velocity (the particle's velocity)
    - initial_position (the particle's initial position)
    - velocity (the particle's velocity)
    - mass (the particle's mass)
    - trajectory (a list describing points in the particle's path)
    """
    def __init__(self, identifier: Any, start_time: float, position: np.array, velocity: np.array, *args, **kwargs):
        """
        The constructor for Particle.

        :param identifier: The unique identifier the particle should have.
        :param start_time: The time the particle starts being simulated at.
        :param position: The phase space position (n-D position and velocity) the particle starts at.
            A (2*n,) dimensional vector, where n is the number of spatial dimensions.
        """
        self.identifier: int = identifier
        self.lost: bool = False
        self.time: float = start_time

        assert len(position) == len(velocity),\
            "Position and velocity must have same length but have shapes %s and %s" % (position.shape, velocity.shape)
        self.number_of_dimensions = len(position)
        self.phase_space_position: np.array = np.concatenate((position, velocity))  # this copies implicitly
        self._initial_tracked_properties: np.array = self.as_array('tracked')
        self.trajectory: List[np.array] = []

        self.detector_hits: Dict[int, List['ParticleDetectorHit']] = {}  # TODO naming is nonagnostic about detectors...
        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            if value != self.number_of_dimensions:
                raise ValueError(f"Number of dimensions changed to {value}, but {self} had {self.number_of_dimensions}")

    @property
    def position(self):
        return self.phase_space_position[:self.number_of_dimensions]

    @position.setter
    def position(self, value):
        self.phase_space_position[:self.number_of_dimensions] = value

    @property
    def velocity(self):
        return self.phase_space_position[self.number_of_dimensions:]

    @velocity.setter
    def velocity(self, value):
        self.phase_space_position[self.number_of_dimensions:] = value

    @property
    def initial_tracked_properties(self):
        return self._initial_tracked_properties

    @cached_property
    @abstractmethod
    def mass(self) -> float:
        pass

    @property
    @abstractmethod
    def tracked_properties(self):
        return [
            ('position', (np.float64, self.number_of_dimensions)),
            ('velocity', (np.float64, self.number_of_dimensions)),
            ('time', np.float64)
        ]

    @property
    @abstractmethod
    def constant_properties(self):
        return [
            ('identifier', np.int32),
            ('mass', np.float64)
        ]

    def _get_props(self, which):
        if which == 'all':
            return self.tracked_properties + self.constant_properties
        elif which == 'tracked':
            return self.tracked_properties
        elif which == 'constant':
            return self.constant_properties

    def as_dtype(self, which) -> np.array:
        """
        Returns the dtype that the array returned by as_array() will have, when called with the same parameters.

        :param which: The part of this particle's state to return. Refer to the docstring about as_array.
        """
        return np.dtype(self._get_props(which))

    def as_array(self, which) -> np.array:
        """
        Returns a NumPy array representation of the current state of this particle. Either the tracked properties,
        or the constant properties, or all (tracked + constant) can be returned as a single NumPy array with a
        structured NumPy datatype.

        :param which: The part of this particle's state to return as an array.
           Must be one of ['tracked', 'constant', 'all'].
             - 'tracked' will refer to the implementation of tracked_properties,
             - 'constant' will refer to the implementation of constant_properties,
             - 'all' will refer to both and concatenate them in order (tracked + constant).
        """
        assert which in ['tracked', 'constant', 'all'],\
            f"which was '{which}' but must be 'tracked', 'constant', or 'all'"
        props = self._get_props(which)
        return np.array([
            tuple(getattr(self, prop[0]) for prop in props)
        ], dtype=self.as_dtype(which))

    @property
    def reached_any_detector(self) -> bool:
        """
        Whether this particle has ever reached any detector.
        """
        return not not self.detector_hits  # `not not` to convert to boolean, to not leak data here

    def __str__(self):
        return f"<{self.__class__.__name__} #{self.identifier}>"

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
