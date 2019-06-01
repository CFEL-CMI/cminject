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
    Particle, Field, Device, Calculator
from cminject.experiment_refactored.definitions.basic import \
    SimpleZBoundary, SimpleZDetector, GaussianSphericalSource, \
    ThermallyConductiveSphericalParticle, CuboidBoundary, empty_interval, TrajectoryCalculator
from cminject.experiment_refactored.fields.fluid_flow_field import FluidFlowField
from cminject.experiment_refactored.visualization.plot import plot_particles


class GravityForceField(Field):
    def __init__(self):
        super().__init__()
        self.g = 9.81  # m/s^2
        self.a = np.array([0, 0, -self.g])  # this force points downwards on the Z axis

    def get_z_boundary(self) -> Tuple[float, float]:
        return empty_interval

    def calculate_acceleration(self,
                               position_and_velocity: np.array,
                               time: float,
                               particle: Particle = None) -> np.array:
        return self.a


class ExampleDevice(Device):
    def __init__(self):
        super().__init__(
            field=GravityForceField(),
            boundary=CuboidBoundary(
                (-0.1, 0.1),
                (-0.1, 0.1),
                (-0.1, 0.1)
            )
        )


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


class ParticleTemperatureCalculator(Calculator):
    def __init__(self, field: FluidFlowField, dt: float):
        self.field: FluidFlowField = field
        self.dt: float = dt

    def calculate(self, particle: ThermallyConductiveSphericalParticle, time: float) -> None:
        if self.field.is_particle_inside(particle):
            t = self.field.calc_particle_temperature(particle, time)
            n_c, t_c = self.field.calc_particle_collisions(particle, self.dt, time)

            particle.collision_temperature = t_c
            particle.temperature = t

            if t <= 77.0 and not np.isfinite(particle.time_to_liquid_n):
                particle.time_to_liquid_n = time
            if t_c <= 77.0 and not np.isfinite(particle.collision_time_to_liquid_n):
                particle.collision_time_to_liquid_n = time


def run_example_experiment(vz, nof_particles, flow_field_filename,
                           track_trajectories=False, do_profiling=False, single_threaded=False):
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
        SimpleZDetector(identifier=0, z_position=-0.052)
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

    calculators = [TrajectoryCalculator()] if track_trajectories else []
    calculators += [
        ParticleTemperatureCalculator(
            field=devices[0].field,
            dt=dt
        )
    ]

    experiment = Experiment(
        devices=devices,
        detectors=detectors,
        sources=sources,
        calculators=calculators,
        t_start=t_start,
        t_end=t_end,
        dt=dt,
        delta_z_end=0.001
    )

    print(f"The total Z boundary of the experiment is {experiment.z_boundary}.")
    print("Running experiment...")
    result = experiment.run(do_profiling=do_profiling, single_threaded=single_threaded)
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
    parser.add_argument('-s', help='Run single threaded? CAUTION: Very slow, only for debugging purposes',
                        action='store_true')
    parser.add_argument('-o', help='Output filename for phase space (hdf5 format)', type=str, required=True)
    parser.add_argument('-f', help='Flow field filename (hdf5 format)', type=str, required=True)
    args = parser.parse_args()

    result_list = run_example_experiment(vz=args.v, nof_particles=args.n,
                                         track_trajectories=args.t, do_profiling=args.p,
                                         single_threaded=args.s, flow_field_filename=args.f)

    # Construct final phase space
    phase_space = [particle.get_phase() for particle in result_list]
    # Construct (detector id -> list[particle phase]) map from all detector hits
    detector_hits = {}
    for particle in result_list:
        if particle.detector_hits:
            for detector_id, hits in particle.detector_hits.items():
                if detector_id not in detector_hits:
                    detector_hits[detector_id] = []
                detector_hits[detector_id] += hits

    # Write final phase space, trajectories, and hits with phases at each detector
    with h5py.File(args.o, 'w') as f:
        f['particles'] = np.array(phase_space)
        # TODO this assumes all particles are the same
        f['particles'].attrs['description'] = result_list[0].get_phase_description()

        for particle in result_list:
            if particle.trajectory:
                f[f'trajectories/{particle.identifier}'] = np.array(particle.trajectory)
                f[f'trajectories/{particle.identifier}'].attrs['description'] = [
                    't', 'x', 'y', 'z', 'u', 'v', 'w'
                ]

        for detector_id, hits in detector_hits.items():
            f[f'detector_hits/{detector_id}'] = np.array([hit.particle_phase for hit in hits])
            f[f'detector_hits/{detector_id}'].attrs['description'] = hits[0].particle_phase_description

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