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
from typing import List, Tuple

import numpy as np
import h5py

from cminject.experiment_refactored.definitions.experiment import Experiment
from cminject.experiment_refactored.definitions.base_classes import \
    Particle, Field, Device
from cminject.experiment_refactored.definitions.basic import \
    infinite_interval, SimpleZBoundary, plot_particles, SimpleZDetector, GaussianSphericalSource,\
    ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.fields.fluid_flow_field import FluidFlowField


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
        self.field: FluidFlowField = FluidFlowField(
            filename=filename,
            density=density, dynamic_viscosity=dynamic_viscosity, scale_slip=scale_slip
        )
        self.boundary: SimpleZBoundary = SimpleZBoundary(self.field.min_z, self.field.max_z)
        super().__init__(field=self.field, boundary=self.boundary)

    def is_particle_inside(self, particle: Particle) -> bool:
        return self.field.is_particle_inside(particle)


def calculate_temperatures(particle: ThermallyConductiveSphericalParticle,
                           devices: List[Device],
                           dt: float,
                           time: float):
    # TODO make this O(1) instead of O(n) by declaring data dependencies explicitly? e.g. a device index?
    for device in devices:
        if isinstance(device, FluidFlowFieldDevice) and device.is_particle_inside(particle):
            t = device.field.calc_particle_temperature(particle, time)
            n_c, t_c = device.field.calc_particle_collisions(particle, dt, time)
            particle.collision_temperature = t_c
            particle.temperature = t

            if t <= 77.0 and particle.time_to_liquid_n is None:
                particle.time_to_liquid_n = time
            if t_c <= 77.0 and particle.collision_time_to_liquid_n is None:
                particle.collision_time_to_liquid_n = time

            # Here, we only care about the first found fluid flow device that a particle is in, so break.
            break


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
    detectors = [
        SimpleZDetector(identifier=0, z_position=-0.052)
    ]

    print(f"The initial velocity of the particles is around {vz} m/s in Z direction.")
    sources = [GaussianSphericalSource(
        nof_particles,
        position=([0.0, 0.0, 0.0048], [0.0002, 0.0002, 0.00001]),
        velocity=([0.0, 0.7, vz], [0.10, 0.30, 2.00]),
        radius=(2.45e-7, 7.5e-9),
        seed=1000,
        rho=1050.0,

        subclass=ThermallyConductiveSphericalParticle,
    )]

    experiment = Experiment(
        devices=devices,
        detectors=detectors,
        sources=sources,
        t_start=0.0,
        t_end=1.0,
        dt=1.0e-6,
        track_trajectories=track_trajectories,
        delta_z_end=0.001,
        calc_additional=calculate_temperatures
    )

    print(f"The total Z boundary of the experiment is {experiment.z_boundary}.")
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
    parser.add_argument('-o', help='Output filename for phase space (hdf5 format)', required=True)
    parser.add_argument('-f', help='Flow field filename (hdf5 format)', type=str, required=True)
    args = parser.parse_args()

    result_list = run_example_experiment(vz=args.v, nof_particles=args.n,
                                         track_trajectories=args.t, do_profiling=args.p,
                                         flow_field_filename=args.f)

    # Construct phase space
    phase_space = []
    trajectories = {}
    for particle in result_list:
        detected_position = None
        if particle.detector_hits:
            detected_position = next(iter(particle.detector_hits.values()))

        if detected_position is not None:
            trajectories[particle.identifier] = np.array(particle.trajectory)
            phase = np.concatenate([
                [particle.radius, particle.mass],
                particle.initial_position,
                particle.initial_velocity,
                detected_position[0][0:2],
                particle.velocity,
                [
                    particle.temperature,
                    particle.collision_temperature,
                    particle.time_to_liquid_n,
                    particle.collision_time_to_liquid_n,
                    particle.time_of_flight,
                    particle.identifier
                ]
            ])
            phase_space.append(phase)

    with h5py.File(args.o, 'w') as f:
        f['particles'] = np.array(phase_space)
        for id, traj in trajectories.items():
            f[f'trajectories/{id}'] = np.array(traj)

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