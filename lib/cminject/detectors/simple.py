#!/usr/bin/env python3
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMInject
#
# Copyright (C) 2018,2020 CFEL Controlled Molecule Imaging group
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any 
# later version. 
#
# If you use this program for scientific work, you should correctly reference it; see the LICENSE.md file for details. 
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

from typing import Optional, Dict, Any, Tuple

import numpy as np

from cminject.base import Detector
from cminject.utils.global_config import ConfigSubscriber, GlobalConfig, ConfigKey


class SimpleZDetector(Detector, ConfigSubscriber):
    def _hit_position(self, phase_space_position: np.array) -> Optional[np.array]:
        raise Exception("This should never be called, the implementation should be swapped out by the constructor!")

    def __init__(self, z_position: float, identifier: str = None):
        """
        Constructs a new SimpleZDetector.
        :param z_position: The z position to put the detector at.
        :param identifier: The identifier of the detector. If None, ("SimpleZ@" + str(z_position)) is used.
        """
        self.z_position = z_position
        self._z_boundary = (z_position, z_position)
        self.particle_distances: Dict[Any, float] = {}

        self._hit_position = None
        self._number_of_dimensions = None
        super().__init__(identifier=identifier if identifier is not None else f'SimpleZ@{z_position}')

        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            number_of_dimensions = value
            self._number_of_dimensions = number_of_dimensions

            if number_of_dimensions == 3:
                self._hit_position = self._hit_position3
            elif number_of_dimensions == 2:
                self._hit_position = self._hit_position2
            else:
                raise ValueError(f"{self.__class__.__name__} requires 2 or 3 spatial dimensions.")

    def _has_particle_reached_detector(self, identifier: str, position: np.array) -> bool:
        if position[-1] == self.z_position:
            return True
        prev_distance = self.particle_distances.get(identifier, None)
        curr_distance = position[-1] - self.z_position
        self.particle_distances[identifier] = curr_distance

        return prev_distance is not None and np.sign(curr_distance) != np.sign(prev_distance)

    def _hit_position3(self, ps_position: np.array) -> Optional[np.array]:
        if ps_position[2] == self.z_position:
            return ps_position[:3]

        x, y, z, u, v, w = ps_position
        d = abs(z) - abs(self.z_position)
        t = d / abs(w)
        x = x - u * t
        y = y - v * t
        return np.array([x, y, self.z_position])

    def _hit_position2(self, ps_position: np.array) -> Optional[np.array]:
        if ps_position[1] == self.z_position:
            return ps_position[:2]

        x, z, vx, vz = ps_position
        d = abs(z) - abs(self.z_position)
        t = d / abs(vz)
        x = x - vx * t
        return np.array([x, self.z_position])

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary
