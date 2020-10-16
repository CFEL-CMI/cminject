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
A collection of subclasses of :class:`cminject.base.Action`. Contains (at least):
  * actions that model Brownian motion by applying a random force after every integration step
  * actions that model other physical changes required by devices (e.g., .temperature.MolecularFlowUpdateTemperature)
  * a simple action to track the particle's trajectory
  * ...your own and more?
"""

from cminject.base import Action, Particle


class TrackTrajectory(Action):
    """
    Appends the array of tracked properties to the trajectory of the particle.
    """
    def __call__(self, particle: Particle, time: float) -> bool:
        # TODO is list appending and converting once to np.array more efficient in the end than appending to np.array?
        particle.trajectory.append(particle.as_array('tracked'))
        return False


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End: