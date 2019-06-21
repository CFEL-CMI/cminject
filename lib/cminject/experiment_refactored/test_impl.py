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
import os
from typing import Tuple

import numpy as np
import h5py
from cminject.experiment_refactored.definitions.result_storage import HDF5ResultStorage

from cminject.experiment_refactored.experiment import Experiment
from cminject.experiment_refactored.definitions.base import \
    Particle, Device
from cminject.experiment_refactored.definitions.property_updaters import TrajectoryPropertyUpdater, \
    ParticleTemperaturePropertyUpdater
from cminject.experiment_refactored.definitions.sources import GaussianSphericalSource
from cminject.experiment_refactored.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.definitions.detectors import SimpleZDetector
from cminject.experiment_refactored.definitions.boundaries import SimpleZBoundary
from cminject.experiment_refactored.fields.fluid_flow_field import StokesFluidFlowField
from cminject.experiment_refactored.visualization.plot import get_3d_figure_from_particles, get_hist2d_figure_from_hdf5


class FluidFlowFieldDevice(Device):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

    def __init__(self, filename: str, density: float, dynamic_viscosity: float, scale_slip: float):
        self.field: StokesFluidFlowField = StokesFluidFlowField(
            filename=filename,
            density=density, dynamic_viscosity=dynamic_viscosity, scale_slip=scale_slip
        )
        self.boundary: SimpleZBoundary = SimpleZBoundary(tuple(self.field._z_boundary))
        super().__init__(field=self.field, boundary=self.boundary)

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.field.is_particle_inside(particle)


def run_example_experiment(vz, nof_particles, flow_field_filename,
                           track_trajectories=False, single_threaded=False):
    print("Setting up experiment...")
    print(f"The initial velocity of the particles is around {vz} m/s in Z direction.")

    t_start, dt, t_end = 0.0, 1.0e-6, 1.0
    print(f"Simulating from t0={t_start} to {t_end} with dt={dt}.")

    devices = [
        FluidFlowFieldDevice(
            filename=flow_field_filename,
            density=0.0009,
            dynamic_viscosity=1.02e-6,
            scale_slip=4.1
        )
    ]
    detectors = [
        SimpleZDetector(identifier=0, z_position=0.0),
        SimpleZDetector(identifier=1, z_position=-0.052)
    ]
    sources = [
        GaussianSphericalSource(
            nof_particles,
            position=([0.0, 0.0, 0.0048], [0.0002, 0.0002, 0.00001]),
            velocity=([0.0, 0.7, vz], [0.10, 0.30, 2.00]),
            radius=(2.45e-7, 7.5e-9),
            seed=1000,
            rho=1050.0,

            subclass=ThermallyConductiveSphericalParticle,
        )
    ]

    property_updaters = [TrajectoryPropertyUpdater()] if track_trajectories else []
    property_updaters += [
        ParticleTemperaturePropertyUpdater(
            field=devices[0].field,
            dt=dt
        )
    ]

    experiment = Experiment(
        devices=devices,
        detectors=detectors,
        sources=sources,
        property_updaters=property_updaters,
        time_interval=(t_start, t_end, dt),
        delta_z_end=0.001
    )

    print(f"The total Z boundary of the experiment is {experiment.z_boundary}.")
    print("Running experiment...")
    experiment_result = experiment.run(single_threaded=single_threaded)
    print("Done running experiment.")
    return experiment_result


def main(args):
    if os.path.isfile(args.output_file):
        raise ValueError("Output file already exists! Please delete or move the existing file.")

    result_particles = run_example_experiment(
        vz=args.z_velocity,
        nof_particles=args.nof_particles,
        flow_field_filename=args.flow_field,
        track_trajectories=args.store_traj,
        single_threaded=args.single_threaded
    )

    HDF5ResultStorage(args.output_file).store_results(result_particles)

    if args.visualize:
        print(f"Plotting {len(result_particles)} particles...")
        fig = get_hist2d_figure_from_hdf5(
            args.output_file,
            [
                ('x', 'y'), ('y', 'z'), ('x', 'z'),
                ('vx', 'vy'), ('vy', 'vz'), ('vx', 'vz')
            ]
        )
        fig.show()
        return fig


if __name__ == '__main__':
    def natural_number(x):
        x = int(x)
        if x < 1:
            raise argparse.ArgumentTypeError("Minimum is 1")
        return x


    parser = argparse.ArgumentParser(prog='cminject',
                                     formatter_class=argparse.MetavarTypeHelpFormatter)
    parser.add_argument('-n', '--nof-particles', help='Number of particles',
                        type=natural_number, required=True)
    parser.add_argument('-v', '--z-velocity', help='Initial average velocity (in Z direction)',
                        type=float, required=True)
    parser.add_argument('-o', '--output-file', help='Output filename for phase space (hdf5 format)',
                        type=str, required=True)
    parser.add_argument('-f', '--flow-field', help='Flow field filename (hdf5 format)',
                        type=str, required=True)
    parser.add_argument('-t', '--store-traj', help='Store trajectories?',
                        action='store_true')
    parser.add_argument('-V', '--visualize', help='Visualize the particles after running?',
                        action='store_true')
    parser.add_argument('-s', '--single-threaded',
                        help='Run single threaded? CAUTION: Very slow, only for debugging purposes',
                        action='store_true')
    main(parser.parse_args())


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""