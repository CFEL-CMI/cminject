#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120; truncate-lines: t -*-
#
# This file is part of CMIstark
# Copyright (C) 2008 Jochen Küpper <jochen.kuepper@cfel.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
"""Unit-tests of Stark effect calculations

Copyright (C) 2008,2014 Jochen Küpper <jochen.kuepper@cfel.de>"""

__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"

import numpy as num
import os
import unittest

import cmiext
import cmiext.convert as convert
from cmiext.state import State

import cmistark.molecule as molecule
import cmistark.starkeffect as starkeffect


class StarkCalculationBenzonitrile(unittest.TestCase):
    """Test the results of Stark effect calculations using the molecular parameters of benzonitrile"""

    @classmethod
    def setUpClass(self):
        """Run before the first test"""
        # set molecular parameters
        self.param = starkeffect.CalculationParameter
        self.param.isomer = 0
        self.param.watson = 'A'
        self.param.symmetry = 'C2a'
        self.param.rotcon = convert.Hz2J(num.array([5655.2654e6, 1546.875864e6, 1214.40399e6]))
        self.param.quartic = convert.Hz2J(num.array([45.6, 938.1, 500, 10.95, 628]))
        self.param.dipole = convert.D2Cm(num.array([4.5152, 0., 0.]))
        # calculation details
        self.param.M = [0, 1]
        self.param.Jmin = 0
        self.param.Jmax_calc = 15
        self.param.Jmax_save =  3
        self.param.dcfields = convert.kV_cm2V_m(num.linspace(0., 100., 5))
        # create Molecule object and specify storage file
        self.param.name = "__cmiext_test_starkeffect"
        self.storagename = self.param.name + ".molecule"
        if os.path.exists(self.storagename):
            raise EnvironmentError("Test storage file already exists, not overwriting")
        self.bn = molecule.Molecule(storage=self.storagename, name=self.param.name)
        # calculate Stark energies
        self.bn.starkeffect_calculation(self.param)

    @classmethod
    def tearDownClass(self):
        """Run after the last test"""
        del(self.bn)
        os.remove(self.storagename)


    def setUp(self):
        """Run before every single test method"""
        pass

    def tearDown(self):
        """Run after every single test method"""
        pass

    def test_fieldfree(self):
        """Test field-free energies (at 0 kV/cm)"""
        # comparing to 0 is dangerous as assertAlmostEqual compares a specified number of digits after 0 i.e. 7
        # and a typical energy is 1e-23 J therefore  convert to Hz before test
        self.assertAlmostEqual(0., convert.J2Hz(self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][0]), 7,
                               "Field-free ground state energy is wrong: expected %g MHz, got %g MHz" \
                               % (convert.J2MHz(0), convert.J2MHz(self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][0])))
        # 101 state has an energy of B+C + centrifugal distortion
        self.assertAlmostEqual(2761.2796716e6,
                               convert.J2Hz(self.bn.starkeffect(State(1, 0, 1, 0, 0))[1][0]), 7,
                               "Field-free ground state energy is wrong: expected %g MHz, got %g MHz" \
                               % (2761.2796716e6, convert.J2MHz(self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][0])))

    def test_hundred(self):
        """Test some state energies at 100 kV/cm

        With our setup, these are the fifth values in the list of fields/energies.
        """
        # test (once) that the fields are correct
        self.assertAlmostEqual(convert.kV_cm2V_m(100.), self.bn.starkeffect(State(0, 0, 0, 0, 0))[0][4], 7,
                               "Field-strength is wrong")
        # test energies for different states at 100 kV/cm
        self.assertAlmostEqual(1., -1.34489847e-22 / self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][4], 7,
                               "Ground state energy is wrong: expected %g MHz, got %g MHz" \
                               % (convert.J2MHz(-1.34489847e-22), convert.J2MHz(self.bn.starkeffect(State(0, 0, 0, 0, 0))[1][4])))


if __name__ == '__main__':
    unittest.main()
