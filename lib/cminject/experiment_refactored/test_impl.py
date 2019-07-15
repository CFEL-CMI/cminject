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
from typing import Tuple, List

from cminject.experiment_refactored.definitions.base import Device
from cminject.experiment_refactored.definitions.boundaries import GridFieldBasedBoundary
from cminject.experiment_refactored.definitions.detectors import SimpleZDetector
from cminject.experiment_refactored.definitions.fields.stokes_fluid_flow_field import StokesFluidFlowField
from cminject.experiment_refactored.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.definitions.property_updaters import TrajectoryPropertyUpdater, \
    MirrorSymmetryPropertyUpdater, BrownianMotionPropertyUpdater
from cminject.experiment_refactored.definitions.result_storage import HDF5ResultStorage
from cminject.experiment_refactored.definitions.sources import VariableDistributionSource
from cminject.experiment_refactored.experiment import Experiment


class StokesFluidFlowFieldDevice(Device):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

    def __init__(self, *args, **kwargs):
        field = StokesFluidFlowField(*args, **kwargs)
        self.fields: List[StokesFluidFlowField] = [field]
        self.boundary: GridFieldBasedBoundary = GridFieldBasedBoundary(field)
        super().__init__(fields=self.fields, boundary=self.boundary)


def run_example_experiment(args):
    if os.path.isfile(args.output_file):
        raise ValueError("Output file already exists! Please delete or move the existing file.")

    print("Setting up experiment...")
    t_start, t_end, dt = args.time_interval

    print(f"Simulating in {args.dimensions}D space from t0={t_start} to {t_end} with dt={dt}.")

    devices = [
        StokesFluidFlowFieldDevice(
            filename=args.flow_field,
            temperature=args.flow_temperature,
            density=args.density,
            dynamic_viscosity=args.flow_dynamic_viscosity,
            m_gas=args.flow_gas_mass,
            slip_correction_scale=args.flow_scale_slip,
        )
    ]

    detectors = [
        SimpleZDetector(identifier=i, z_position=pos) for i, pos in enumerate(args.detectors)
    ]

    sources = [
        VariableDistributionSource(
            args.nof_particles,
            position=args.position,
            velocity=args.velocity,
            radius=args.radius,
            rho=args.density,

            seed=args.seed,
            subclass=ThermallyConductiveSphericalParticle
        ),
    ]
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
        property_updaters += [BrownianMotionPropertyUpdater(field=devices[0].fields[0], dt=dt)]
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
    HDF5ResultStorage(args.output_file).store_results(result_particles)
    print(f"Saved results to {args.output_file}.")

    return result_particles


# Everything below is just argument parsing even if it looks like a lot of code :)
if __name__ == '__main__':
    def dist_description(x):
        values = x.split()
        kind = values[0]
        if kind == 'G':
            if len(values) != 3:
                parser.error("Need a mu and sigma for a gaussian distribution description!")
            return {'kind': 'gaussian', 'mu': float(values[1]), 'sigma': float(values[2])}
        elif kind == 'L':
            if len(values) != 3:
                parser.error("Need a min and max for a linear distribution description!")
            return {'kind': 'linear', 'min': float(values[1]), 'max': float(values[2])}
        else:
            parser.error(f"Unknown distribution kind: {kind}!")

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

    parser.add_argument('-r', '--radius', help='Distribution description for the radius.',
                        type=dist_description)
    parser.add_argument('-p', '--position',
                        help='Distribution description for the position.',
                        nargs='*', type=dist_description)
    parser.add_argument('-v', '--velocity',
                        help='Distribution description for the velocity.',
                        nargs='*', type=dist_description)
    parser.add_argument('-rho', '--density', help='Density of the particles.',
                        type=float)

    parser.add_argument('-fG', '--flow-gas', choices=['N', 'He'],
                        help='Shortcut for -fv, -fd and -fg for common gas types. Warns if values are unknown'
                             'for the given flow temperature.',
                        type=str)
    parser.add_argument('-ft', '--flow-temperature', help='Temperature of the gas in the flow field',
                        type=float)
    parser.add_argument('-fv', '--flow-dynamic-viscosity', help='Dynamic viscosity of the flow field',
                        type=float)
    parser.add_argument('-fd', '--flow-density', help='Density of the flow field',
                        type=float)
    parser.add_argument('-fs', '--flow-scale-slip', help='Slip scale factor of the flow field',
                        type=float)
    parser.add_argument('-fm', '--flow-gas-mass', help='Mass of the gas particles in the flow field',
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

    parser.set_defaults(dimensions=3,
                        time_interval=[0.0, 1.8, 1e-6],
                        seed=1000,
                        # Assume Helium at 4.0 K by default
                        flow_dynamic_viscosity=1.02e-6, flow_temperature=4.0, flow_gas_mass=6.6e-27,
                        flow_density=0.0009, flow_scale_slip=4.1,
                        radius={'kind': 'gaussian', 'mu': 2.45e-7, 'sigma': 7.5e-9},
                        position=[
                            {'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.0002},
                            {'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.0002},
                            {'kind': 'gaussian', 'mu': 0.0048, 'sigma': 0.00001}
                        ],
                        velocity=[
                            {'kind': 'gaussian', 'mu': 0.0, 'sigma': 0.1},
                            {'kind': 'gaussian', 'mu': 0.7, 'sigma': 0.3},
                            {'kind': 'gaussian', 'mu': -43.0, 'sigma': 2.0},
                        ],
                        density=1050.0,
                        detectors=[0.0, -0.052])
    args = parser.parse_args()

    # Set mu/m_gas if --flow-gas was set
    if args.flow_gas is not None:
        print(f"Automatically setting mu and m_gas for gas type '{args.flow_gas}' at {args.flow_temperature}K.")
    if args.flow_gas == 'He':
        args.flow_gas_mass = 6.6e-27
        if args.flow_temperature == 4.0:
            args.flow_dynamic_viscosity = 1.02e-6
        else:
            args.flow_dynamic_viscosity = 1.96e-5
            if args.flow_temperature != 293.15:
                print(f"WARNING: Unknown mu/m for gas type {args.flow_gas} at {args.flow_temperature}. Using "
                      f"approximation for room temperature (293.15K).")
    elif args.flow_gas == 'N':
        args.flow_gas_mass = 4.7e-26
        args.flow_dynamic_viscosity = 1.76e-5
        if args.flow_temperature != 293.15:
            print(f"WARNING: Unknown mu/m for gas type {args.flow_gas} at {args.flow_temperature}. Using "
                  f"approximation for room temperature (293.15K).")

    # Verify dimensionality match for position description and dimensions parameter
    if len(args.position) != args.dimensions:
        parser.error("The length of the position description vectors must match the simulation dimensionality!")
    if len(args.velocity) != args.dimensions:
        parser.error("The length of the velocity description vectors must match the simulation dimensionality!")

    run_example_experiment(args)


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
