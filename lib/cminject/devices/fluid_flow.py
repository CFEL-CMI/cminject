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

"""A device for one flow field."""


from enum import Enum

from cminject.actions.brownian_motion import MolecularFlowBrownianMotionStep, StokesBrownianMotionStep
from cminject.actions.temperature import MolecularFlowUpdateTemperature
from cminject.base import Device
from cminject.boundaries.field_based import GridFieldBasedBoundary
from cminject.fields.fluid_flow import MolecularFlowDragForceField, StokesDragForceField


class FlowType(Enum):
    """An enumeration representing the available types of fluid flow models."""
    MOLECULAR_FLOW = 'Molecular flow (modified Epstein force by Nils Roth)'
    STOKES = 'Stokes'


class FluidFlowDevice(Device):
    """
    A device for one flow field. It derives its boundary automatically based on the minimal cuboid enclosing
    the flow field.

    For calculating the force, there are two options available (see the FlowType definition): The Stokes drag force
    model (for Kn << 1), and a modified Epstein force (for Kn >> 1), which was developed by Nils Roth.
    """
    def __init__(self, filename: str, flow_type: FlowType, brownian_motion: bool, *args, **kwargs):
        field_class = StokesDragForceField if flow_type == FlowType.STOKES else MolecularFlowDragForceField
        bm_class = StokesBrownianMotionStep if flow_type == FlowType.STOKES else MolecularFlowBrownianMotionStep

        field = field_class(filename=filename, *args, **kwargs)
        boundary = GridFieldBasedBoundary(field=field)
        actions = []
        if flow_type == FlowType.MOLECULAR_FLOW:
            actions.append(MolecularFlowUpdateTemperature(field))
        if brownian_motion:
            actions.append(bm_class(field))

        super().__init__(fields=[field], boundary=boundary, actions=actions)
