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

import argparse
from typing import List, Tuple, Optional

import numpy as np

from cminject.experiment_refactored.experiment import Experiment
from cminject.experiment_refactored.base_classes import \
    Particle, Field, Device, ZBoundedMixin
from cminject.experiment_refactored.basic import \
    infinite_interval, SimpleZBoundary, plot_particles, SimpleZDetector, GaussianSphericalSource, InfiniteBoundary, \
    CuboidBoundary
from cminject.experiment_refactored.fluid_flow_field import FluidFlowField


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


class FluidFlowFieldDevice(Device):
    def __init__(self, filename: str, density: float, dynamic_viscosity: float, scale_slip: float):
        super().__init__(
            field=FluidFlowField(
                filename=filename,
                density=density, dynamic_viscosity=dynamic_viscosity, scale_slip=scale_slip
            ),
            boundary=CuboidBoundary(
                (-0.01, 0.01),
                (-0.01, 0.01),
                (-0.06, 0.02),
            )
        )

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.boundary.is_particle_inside(particle) or self.field.is_particle_inside(particle)


def run_example_experiment(vz, nof_particles, flow_field_filename, track_trajectories=False, do_profiling=False):
    print("Setting up experiment...")
    devices = [
        # ExampleDevice(),
        FluidFlowFieldDevice(
            filename=flow_field_filename,
            density=0.0009,
            dynamic_viscosity=1.02e-6,
            scale_slip=4.1
        )
    ]
    detectors = [SimpleZDetector(identifier=0, z_position=-0.052)]
    sources = [GaussianSphericalSource(
        nof_particles,
        position=([0.0, 0.0, 0.0048], [0.0002, 0.0002, 0.00001]),
        velocity=([0.0, 0.7, vz], [0.10, 0.30, 2.00]),
        radius=(2.45e-7, 7.5e-9),
        seed=1000,
        rho=1050.0
    )]
    experiment = Experiment(
        devices=devices,
        detectors=detectors,
        sources=sources,
        t_start=0.0,
        t_end=1.0,
        dt=1.0e-6,
        track_trajectories=track_trajectories
    )

    print("Running experiment...")
    result = experiment.run(do_profiling=do_profiling, single_threaded=False)
    print("Done running experiment.")
    return result


def natural_number(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("Minimum is 1")
    return x


def main():
    parser = argparse.ArgumentParser(prog='cminject',
                                     formatter_class=argparse.MetavarTypeHelpFormatter)
    parser.add_argument('-n', help='Number of particles', type=natural_number, required=True)
    parser.add_argument('-v', help='Initial average velocity (in Z direction)', type=float, required=True)
    parser.add_argument('-t', help='Store trajectories?', action='store_true')
    parser.add_argument('-p', help='Do profiling? (CAUTION: generates lots of profiling dump files)',
                        action='store_true')
    parser.add_argument('-f', help='Flow field filename (txt format)', type=str, required=True)
    args = parser.parse_args()

    result_list = run_example_experiment(vz=-10.0, nof_particles=args.n,
                                         track_trajectories=args.t, do_profiling=args.p,
                                         flow_field_filename=args.f)
    print(f"Plotting {len(result_list)} particles...")
    plot_particles(result_list, plot_trajectories=args.t)


if __name__ == '__main__':
    main()



"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""