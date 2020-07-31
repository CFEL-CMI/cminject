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

from typing import Tuple, Union

from cminject.base import Device
from cminject.boundaries.simple import CuboidBoundary
from cminject.fields.fluid_flow import StokesDragForceField
from cminject.fields.photophoresis import DesyatnikovPhotophoreticLaserField


class DesyatnikovPhotophoresisDevice(Device):
    """
    A Device for simulating photophoretic effects with the model described in: A. Desyatnikov, V. Shvedov, A. Rode,
    W. Krolikowski, and Y. Kivshar, "Photophoretic manipulation of absorbing aerosol particles with vortex beams: theory
    versus experiment," Opt. Express  17, 8201-8211 (2009).
    """
    def __init__(self, r_boundary: Tuple[float, float], z_boundary: Tuple[float, float],
                 gas_temperature: float, gas_viscosity: float, gas_thermal_conductivity: float,
                 gas_density: Union[float, Tuple[float, StokesDragForceField]], gas_mass: float,
                 beam_power: float, beam_waist_radius: float,  z_position: float = 0.0):
        """
        The constructor for DesyatnikovPhotophoresisDevice.

        :param r_boundary: The r interval where the laser is active
        :param z_boundary: The z interval where the laser is active
        :param gas_temperature: The temperature of the gas the laser is positioned in [K]
        :param gas_viscosity: The viscosity of the gas the laser is positioned in [Pa*s]
        :param gas_thermal_conductivity: The thermal conductivity of the gas the laser is positioned in [W/(m*K)]
        :param gas_density: The density of the gas the laser is positioned in, either as a fixed value or a tuple of a
            default value (a float) and a flow field to base the density on an interpolated pressure.
        :param gas_mass: The mass of a single gas molecule [kg]
        :param beam_power: The total power of the laser [W]
        :param beam_waist_radius: The radius of the laser beam waist [m]
        :param z_position: The Z position of the laser beam waist (0.0 by default)
        """
        pp_field = DesyatnikovPhotophoreticLaserField(
            gas_viscosity=gas_viscosity, gas_temperature=gas_temperature,
            gas_thermal_conductivity=gas_thermal_conductivity, gas_density=gas_density, gas_mass=gas_mass,
            beam_power=beam_power, beam_waist_radius=beam_waist_radius, z_position=z_position
        )
        boundary: CuboidBoundary = CuboidBoundary(intervals=[r_boundary, z_boundary])
        super().__init__(fields=[pp_field], boundary=boundary)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
