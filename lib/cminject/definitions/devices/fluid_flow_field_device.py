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

from enum import Enum

import numpy as np

from cminject.definitions.devices.base import Device
from cminject.definitions.particles.base import Particle
from cminject.definitions.boundaries.grid_field_based import GridFieldBasedBoundary
from cminject.definitions.fields.fluid_flow import StokesDragForceField, MolecularFlowDragForceField


class FlowType(Enum):
    MOLECULAR_FLOW = 'Molecular flow (modified Epstein force by Nils Roth)'
    STOKES = 'Stokes'


class FluidFlowFieldDevice(Device):
    """
    A device for one flow field. It derives its boundary automatically based on the minimal cuboid enclosing
    the flow field.

    For calculating the force, there are two options available (see the FlowType definition): The Stokes drag force
    model (for Kn << 1), and a modified Epstein force (for Kn >> 1), which was developed by Nils Roth.
    """
    def __init__(self, flow_type: FlowType, filename: str, *args, **kwargs):
        if flow_type is FlowType.MOLECULAR_FLOW:
            field = MolecularFlowDragForceField(filename=filename, *args, **kwargs)
        else:
            field = StokesDragForceField(filename=filename, *args, **kwargs)

        boundary = GridFieldBasedBoundary(field=field)
        super().__init__(fields=[field], boundary=boundary)

    def calculate_acceleration(self, particle: Particle, time: float) -> np.array:
        """
        This method is overridden since we know this device only has one field -- so we can avoid the overhead
        of summing a list of 1 item that the default implementation of calculate_acceleration of Device has
        """
        return self.fields[0].calculate_acceleration(particle, time)
