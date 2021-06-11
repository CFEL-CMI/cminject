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
import os

import numpy as np
import pytest

from cminject.particles import SphericalParticle
from cminject.utils.global_config import GlobalConfig, ConfigKey

__PARTICLE_IDENTIFIER = 0


@pytest.fixture(scope="function")
def spherical_particle_2d():
    """
    Simple fixture to return a sample SphericalParticle instance.
    The radius is 13.5nm, density is 19320kg/m^3 (gold), and position
    and velocity are initialized to both be [0.0, 0.0].

    Increases the returned particle identifier after each usage.
    """
    global __PARTICLE_IDENTIFIER
    particle = SphericalParticle(
        identifier=__PARTICLE_IDENTIFIER,
        start_time=0.0, radius=13.5e-9, rho=19320,
        # choose some coordinate that's not on a grid point
        position=np.array([0.0, 0.0], dtype=np.float64),
        velocity=np.array([0.0, 0.0], dtype=np.float64)
    )
    __PARTICLE_IDENTIFIER += 1
    return particle


@pytest.fixture(scope="function")
def dims2d():
    GlobalConfig().set(ConfigKey.NUMBER_OF_DIMENSIONS, 2)


@pytest.fixture(scope="function")
def dt1e_5():
    GlobalConfig().set(ConfigKey.TIME_STEP, 1e-5)


@pytest.fixture(scope="function")
def example_field_file():
    return os.path.join(
        os.path.dirname(__file__),
        '..', 'examples', 'example_field.h5'
    )
