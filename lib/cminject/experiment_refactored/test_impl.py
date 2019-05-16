"""
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
"""

from typing import List, Tuple, Optional

import numpy as np

from cminject.experiment_refactored.experiment import Experiment
from cminject.experiment_refactored.base_classes import\
    Particle, Field, Device
from cminject.experiment_refactored.basic import\
    infinite_interval, SimpleZBoundary, plot_trajectories, SimpleZDetector, GaussianSphericalSource


class GravityForceField(Field):
    def __init__(self):
        super().__init__()
        self.g = -9.81 * 20  # m/s^2

    def get_z_boundary(self) -> Tuple[float, float]:
        return infinite_interval

    def calculate_acceleration(self,
                               position_and_velocity: np.array,
                               time: float,
                               particle: Particle = None) -> np.array:
        return np.array([
            0, 0, self.g
        ])


class ExampleDevice(Device):
    def __init__(self):
        super().__init__(field=GravityForceField(), boundary=SimpleZBoundary(0.0, 1.0))


def run_example_experiment(vz, do_profiling=False):
    devices = [ExampleDevice()]
    detectors = [SimpleZDetector(identifier=1, z_position=0.5)]
    sources = [GaussianSphericalSource(
        10,
        position=([0.0, 0.0, 0.1], [0.0002, 0.0002, 0.0000]),
        velocity=([0.1, 0.0, vz], [0.0002, 0.0002, 0.00000]),
        radius=(0.3, 0.2),
        seed=100,
        rho=1.0
    )]
    experiment = Experiment(
        devices=devices,
        detectors=detectors,
        sources=sources,
        t_end=0.1,
        dt=0.00001,
        track_trajectories=True
    )
    result = experiment.run(do_profiling=do_profiling)
    return result


def main():
    import sys
    do_profiling = (len(sys.argv) > 1 and sys.argv[1] == 'profile')
    result_list = run_example_experiment(vz=13.0, do_profiling=do_profiling)
    plot_trajectories(result_list)


if __name__ == '__main__':
    main()



"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""