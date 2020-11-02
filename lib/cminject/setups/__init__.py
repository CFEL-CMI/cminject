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

"""
Contains a collection of useful subclasses of :class:`cminject.base.Setup`, to be

  * simply used in simulations, or
  * changed/extended/improved and possibly merged back into this repository via a pull request

New Setup subclasses that are useful to others will gladly be accepted via pull requests as well.
"""

__all__ = ['OneFlowFieldSetup', 'DesyatnikovPhotophoresisSetup']

from .one_flow_field import OneFlowFieldSetup
from .desyatnikov_photophoresis import DesyatnikovPhotophoresisSetup
