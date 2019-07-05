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
    MirrorSymmetryPropertyUpdater, BrownianMotionPropertyUpdater
from cminject.experiment_refactored.definitions.sources import GaussianDistributedSphericalSource, \
    LinearlyDistributedSphericalSource
from cminject.experiment_refactored.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.definitions.detectors import SimpleZDetector
from cminject.experiment_refactored.definitions.boundaries import SimpleZBoundary
from cminject.experiment_refactored.definitions.fields.stokes_fluid_flow_field import StokesFluidFlowField


class StokesFluidFlowFieldDevice(Device):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

    def __init__(self, filename: str, temperature: float,
                 dynamic_viscosity: float, slip_correction_scale: float, m_gas: float):
        self.field: StokesFluidFlowField = StokesFluidFlowField(
            filename=filename,
            slip_correction_scale=slip_correction_scale,
            dynamic_viscosity=dynamic_viscosity,
            m_gas=m_gas,
            temperature=temperature
        )
        self.boundary: SimpleZBoundary = SimpleZBoundary(tuple(self.field.z_boundary))
        super().__init__(field=self.field, boundary=self.boundary)

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.field.is_particle_inside(particle)


def run_example_experiment(args):
    if os.path.isfile(args.output_file):
        raise ValueError("Output file already exists! Please delete or move the existing file.")

    print("Setting up experiment...")
    t_start, t_end, dt = args.time_interval

    print(f"Simulating in {args.dimensions}D space from t0={t_start} to {t_end} with dt={dt}.")
    print(f"The initial velocity of the particles is around {args.velocity[args.dimensions - 1]} m/s in Z direction.")

    devices = [
        StokesFluidFlowFieldDevice(
            filename=args.flow_field,
            temperature=args.flow_temperature,
            dynamic_viscosity=args.flow_dynamic_viscosity,
            m_gas=args.flow_gas_mass,
            slip_correction_scale=args.flow_scale_slip,
        )
    ]

    detectors = [
        SimpleZDetector(identifier=i, z_position=pos) for i, pos in enumerate(args.detectors)
    ]

    sources = [
        GaussianDistributedSphericalSource(
            args.nof_particles,
            position=(args.position, args.position_sigma),
            velocity=(args.velocity, args.velocity_sigma),
            radius=(args.radius, args.radius_sigma),
            seed=args.seed,
            rho=args.density,

            subclass=ThermallyConductiveSphericalParticle,
        )
    ]
    # TODO implement arguments to specify these distributions
    """sources = [
        LinearlyDistributedSphericalSource(
            args.nof_particles,
            position_intervals=[(0, 3e-3), (-0.0253, -0.0253)],
            velocity_intervals=[(0.0, 0.0), (.1, .1)],
            radius_interval=(25e-9, 25e-9),
            rho=args.density,
            seed=args.seed
        )
    ]"""

    property_updaters = []
    if args.brownian:
        property_updaters += [BrownianMotionPropertyUpdater(field=devices[0].field, dt=dt)]
    if args.dimensions == 2:
        property_updaters += [MirrorSymmetryPropertyUpdater()]
    property_updaters += [TrajectoryPropertyUpdater()] if args.store_traj else []

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
    print("Done running experiment. Storing results...")
    print([(p.initial_position, p.position) for p in result_particles])
    HDF5ResultStorage(args.output_file).store_results(result_particles)
    print(f"Saved results to {args.output_file}.")

    return result_particles


if __name__ == '__main__':
    def natural_number(x):
        x = int(x)
        if x < 1:
            raise argparse.ArgumentTypeError("Minimum is 1")
        return x

    parser = argparse.ArgumentParser(prog='cminject',
                                     formatter_class=argparse.MetavarTypeHelpFormatter)
    parser.add_argument('-n', '--nof-particles', help='The number of particles to simulate',
                        type=natural_number, required=True)
    parser.add_argument('-t', '--time-interval',
                        help='The time interval to run the simulation within, as a 3-list (t_start, t_end, dt)',
                        type=float, nargs=3)
    parser.add_argument('-f', '--flow-field', help='Flow field filename (hdf5 format)',
                        type=str, required=True)
    parser.add_argument('-d', '--detectors', help='The Z positions of the detectors', nargs='+',
                        type=float, required=True)
    parser.add_argument('-o', '--output-file', help='Output filename for phase space (hdf5 format)',
                        type=str, required=True)

    parser.add_argument('-r', '--radius', help='Radius mu (Gaussian mean) of the particles',
                        type=float)
    parser.add_argument('-rsg', '--radius_sigma', help='Radius sigma (Gaussian standard deviation) of the particles',
                        type=float)
    parser.add_argument('-p', '--position',
                        help='Initial position mu (Gaussian mean), one for each dimension',
                        nargs='*', type=float)
    parser.add_argument('-psg', '--position_sigma',
                        help='Initial position sigma (Gaussian standard deviation), one for each dimension.',
                        nargs='*', type=float)
    parser.add_argument('-v', '--velocity',
                        help='Initial position mu (Gaussian mean), one for each dimension.',
                        nargs='*', type=float)
    parser.add_argument('-vsg', '--velocity_sigma',
                        help='Initial position sigma (Gaussian standard deviation), one for each dimension.',
                        nargs='*',  type=float)
    parser.add_argument('-rho', '--density', help='Density of the particles',
                        type=float)

    parser.add_argument('-fv', '--flow-dynamic-viscosity', help='Dynamic viscosity of the flow field',
                        type=float)
    parser.add_argument('-fs', '--flow-scale-slip', help='Slip scale factor of the flow field',
                        type=float)
    parser.add_argument('-fm', '--flow-gas-mass', help='Mass of the gas particles in the flow field',
                        type=float)
    parser.add_argument('-ft', '--flow-temperature', help='Temperature of the flow field',
                        type=float)

    parser.add_argument('-T', '--store-traj', help='Store trajectories?',
                        action='store_true')
    parser.add_argument('-s', '--single-threaded',
                        help='Run single threaded? CAUTION: Very slow, only for debugging purposes',
                        action='store_true')
    parser.add_argument('-D', '--dimensions', help='# of spatial dimensions',
                        type=int)
    parser.add_argument('-S', '--seed', help='Seed for the random generator',
                        type=int)
    parser.add_argument('-B', '--brownian', help='Enable Brownian motion', action='store_true')

    parser.set_defaults(density=1050.0, dimensions=3, seed=1000,
                        flow_dynamic_viscosity=1.02e-6, flow_density=0.0009, flow_scale_slip=4.1,
                        radius=2.45e-7, radius_sigma=7.5e-9,
                        position=[0.0, 0.0, 0.0048], position_sigma=[0.0002, 0.0002, 0.00001],
                        velocity=[0.0, 0.7, -43.0], velocity_sigma=[0.10, 0.30, 2.00],
                        detectors=[0.0, -0.052], time_interval=[0.0, 1.0, 1e-6])
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
