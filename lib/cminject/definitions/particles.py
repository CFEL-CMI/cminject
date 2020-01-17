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

from typing import List

import numpy as np
from cminject.definitions.base import Particle


class SphericalParticle(Particle):
    """
    A simple spherical particle that has a radius and a density (rho).
    """
    def __init__(self, *args, radius: float, rho: float, **kwargs):
        super().__init__(*args, **kwargs)
        self.radius = radius
        self.rho = rho
        self.mass = self.rho * 4 / 3 * np.pi * (self.radius ** 3)
        self.number_of_dimensions = None

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions in [1, 2, 3]:
            super().set_number_of_dimensions(number_of_dimensions)
        else:
            raise ValueError("SphericalParticles can only be simulated in 1-, 2-, and 3D space.")

    @property
    def properties(self) -> np.array:
        return np.array([
            self.identifier,
            self.time_of_flight,
            self.radius,
            self.mass
        ])

    @property
    def properties_description(self) -> List[str]:
        return ['ID', 't', 'r', 'm']


class ThermallyConductiveSphericalParticle(SphericalParticle):
    """
    Lacking a better name (so far), this class extends on the SphericalParticle with a thermal conductivity
    that is required by e.g. the photophoretic force and for thermal calculations.
    """
    def __init__(self, *args,
                 thermal_conductivity: float, temperature: float, specific_heat: float,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.thermal_conductivity = thermal_conductivity
        self.temperature = temperature
        self.specific_heat = specific_heat

    @property
    def properties(self) -> np.array:
        return np.concatenate([super().properties, [self.temperature]])

    @property
    def properties_description(self) -> List[str]:
        return super().properties_description + ['T']

class Molecule():
    """
    A molecule that has a mass and a quantum state charchterized by the quantum numbers:
    J, Ka, Kc, m, Isomer
    """
    def __init__(self, mass: float, q_n: dict):

        self.mass = mass
        self.q_n = q_n

    @property
    def properties(self) -> np.array:
        return np.array([
            self.time_of_flight,
            self.mass
        ])

if __name__ == "__main__":
    d = {'J': 2, 'Ka': 2, 'Kc': 2, 'M': 2, 'Isomer':0}

    print('m'+str(d['J']))
    mol = Molecule(34, d)
    print(mol.mass)
    print(mol.q_n)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
