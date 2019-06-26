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

from cminject.experiment_refactored.definitions.result_storage import HDF5ResultStorage

from cminject.experiment_refactored.experiment import Experiment
from cminject.experiment_refactored.definitions.base import \
    Particle, Device
from cminject.experiment_refactored.definitions.property_updaters import TrajectoryPropertyUpdater, \
    ParticleTemperaturePropertyUpdater, RadialSymmetryPropertyUpdater
from cminject.experiment_refactored.definitions.sources import GaussianSphericalSource
from cminject.experiment_refactored.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.definitions.detectors import SimpleZDetector
from cminject.experiment_refactored.definitions.boundaries import SimpleZBoundary
from cminject.experiment_refactored.definitions.fields.stokes_fluid_flow_field import StokesFluidFlowField


class StokesFluidFlowFieldDevice(Device):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

    def __init__(self, filename: str, density: float, dynamic_viscosity: float, scale_slip: float):
        self.field: StokesFluidFlowField = StokesFluidFlowField(
            filename=filename,
            density=density, dynamic_viscosity=dynamic_viscosity, scale_slip=scale_slip
        )
        self.boundary: SimpleZBoundary = SimpleZBoundary(tuple(self.field.z_boundary))
        super().__init__(field=self.field, boundary=self.boundary)

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.field.is_particle_inside(particle)


def run_example_experiment(args):
    if os.path.isfile(args.output_file):
        raise ValueError("Output file already exists! Please delete or move the existing file.")

    print("Setting up experiment...")
    print(f"The initial velocity of the particles is around {args.velocity[args.dimensions - 1]} m/s in Z direction.")

    t_start, dt, t_end = 0.0, 1.0e-6, 1.0
    print(f"Simulating from t0={t_start} to {t_end} with dt={dt}.")

    devices = [
        StokesFluidFlowFieldDevice(
            filename=args.flow_field,
            density=args.flow_density,
            dynamic_viscosity=args.flow_dynamic_viscosity,
            scale_slip=args.flow_scale_slip,
        )
    ]
    detectors = [
        SimpleZDetector(identifier=i, z_position=-0.005 * i)
        for i in range(12)
    ]

    sources = [
        GaussianSphericalSource(
            args.nof_particles,
            position=(args.position, args.position_sigma),
            velocity=(args.velocity, args.velocity_sigma),
            radius=(args.radius, args.radius_sigma),
            seed=args.seed,
            rho=args.density,

            subclass=ThermallyConductiveSphericalParticle,
        )
    ]

    property_updaters = [TrajectoryPropertyUpdater()] if args.store_traj else []
    if 0 == 1 % 2:
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
        delta_z_end=0.001,
        number_of_dimensions=args.dimensions
    )

    print(f"The total Z boundary of the experiment is {experiment.z_boundary}.")
    print("Running experiment...")
    result_particles = experiment.run(single_threaded=args.single_threaded)
    print("Done running experiment.")

    HDF5ResultStorage(args.output_file).store_results(result_particles)
    return result_particles


if __name__ == '__main__':
    def natural_number(x):
        x = int(x)
        if x < 1:
            raise argparse.ArgumentTypeError("Minimum is 1")
        return x

    def float_list(x):
        entries = x.split(',')
        return [float(entry) for entry in entries]

    parser = argparse.ArgumentParser(prog='cminject',
                                     formatter_class=argparse.MetavarTypeHelpFormatter)
    parser.add_argument('-n', '--nof-particles', help='Number of particles',
                        type=natural_number, required=True)
    parser.add_argument('-f', '--flow-field', help='Flow field filename (hdf5 format)',
                        type=str, required=True)
    parser.add_argument('-o', '--output-file', help='Output filename for phase space (hdf5 format)',
                        type=str, required=True)

    parser.add_argument('-r', '--radius', help='Radius mu (Gaussian mean) of the particles',
                        type=float)
    parser.add_argument('-rsg', '--radius_sigma', help='Radius sigma (Gaussian standard deviation) of the particles',
                        type=float)
    parser.add_argument('-d', '--density', help='Density of the particles',
                        type=float)
    parser.add_argument('-p', '--position',
                        help='Initial position mu (Gaussian mean), as a tuple of floats separated by commas.',
                        type=float_list)
    parser.add_argument('-psg', '--position_sigma',
                        help='Initial position sigma (Gaussian standard deviation), as a tuple of floats '
                             'separated by commas.',
                        type=float_list)
    parser.add_argument('-v', '--velocity',
                        help='Initial position mu (Gaussian mean), as a tuple of floats separated by commas.',
                        type=float_list)
    parser.add_argument('-vsg', '--velocity_sigma',
                        help='Initial position sigma (Gaussian standard deviation), '
                             'as a tuple of floats separated by commas. ',
                        type=float_list)

    parser.add_argument('-fd', '--flow-density', help='Density of the gas in the flow field',
                        type=float)
    parser.add_argument('-fv', '--flow-dynamic-viscosity', help='Dynamic viscosity of the flow field',
                        type=float)
    parser.add_argument('-fs', '--flow-scale-slip', help='Slip scale factor of the flow field',
                        type=float)

    parser.add_argument('-t', '--store-traj', help='Store trajectories?',
                        action='store_true')
    parser.add_argument('-s', '--single-threaded',
                        help='Run single threaded? CAUTION: Very slow, only for debugging purposes',
                        action='store_true')
    parser.add_argument('-D', '--dimensions', help='# of spatial dimensions',
                        type=int)
    parser.add_argument('-S', '--seed', help='Seed for the random generator',
                        type=int)

    parser.set_defaults(density=1050.0, dimensions=3, seed=1000,
                        flow_dynamic_viscosity=1.02e-6, flow_density=0.0009, flow_scale_slip=4.1,
                        radius=2.45e-7, radius_sigma=7.5e-9,
                        position=[0.0, 0.0, 0.0048], position_sigma=[0.0002, 0.0002, 0.00001],
                        velocity=[0.0, 0.7, -43.0], velocity_sigma=[0.10, 0.30, 2.00])
    args = parser.parse_args()

    # Verify dimensionality match for position description and dimensions parameter
    if len(args.position) != args.dimensions or len(args.position_sigma) != args.dimensions:
        parser.error("The length of the position description vectors must match the simulation dimensionality!")
    if len(args.velocity) != args.dimensions or len(args.velocity_sigma) != args.dimensions:
        parser.error("The length of the velocity description vectors must match the simulation dimensionality!")

    run_example_experiment(args)


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""