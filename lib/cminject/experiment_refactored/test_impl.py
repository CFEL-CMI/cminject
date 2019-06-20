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
from typing import Tuple

import numpy as np
import h5py

from cminject.experiment_refactored.experiment import Experiment
from cminject.experiment_refactored.definitions.base import \
    Particle, Device
from cminject.experiment_refactored.definitions.property_updaters import TrajectoryPropertyUpdater, \
    ParticleTemperaturePropertyUpdater
from cminject.experiment_refactored.definitions.sources import GaussianSphericalSource
from cminject.experiment_refactored.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.definitions.detectors import SimpleZDetector
from cminject.experiment_refactored.definitions.boundaries import SimpleZBoundary
from cminject.experiment_refactored.fields.fluid_flow_field import FluidFlowField
from cminject.experiment_refactored.visualization.plot import plot_particles


class FluidFlowFieldDevice(Device):
    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self.boundary.z_boundary

    def __init__(self, filename: str, density: float, dynamic_viscosity: float, scale_slip: float):
        self.field: FluidFlowField = FluidFlowField(
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


def main(vz, nof_particles, output_file, flow_field, track_trajectories, single_threaded, visualize=False):
    result_particles = run_example_experiment(vz=vz, nof_particles=nof_particles, flow_field_filename=flow_field,
                                              track_trajectories=track_trajectories, single_threaded=single_threaded)

    # Construct final phase space
    phase_space = [np.concatenate([particle.position, particle.properties]) for particle in result_particles]
    phase_space_description = result_particles[0].position_description + result_particles[0].properties_description

    # Write final phase space, trajectories, and hits with phases at each detector
    with h5py.File(output_file, 'w') as f:
        f['particles'] = np.array(phase_space)
        f['particles'].attrs['description'] = phase_space_description

        for particle in result_particles:
            if particle.trajectory:
                f[f'trajectories/{particle.identifier}'] = np.array(particle.trajectory)
                f[f'trajectories/{particle.identifier}'].attrs['description'] = [
                    't', 'x', 'y', 'z', 'u', 'v', 'w'
                ]

        # Construct inverse association
        # (detector ID -> list of particle hits) from (particle ID -> list of detector hits), as this is what we
        # would like to store in the file
        detector_hits = {}
        for particle in result_particles:
            if particle.detector_hits:
                for detector_id, hits in particle.detector_hits.items():
                    if detector_id not in detector_hits:
                        detector_hits[detector_id] = []
                    detector_hits[detector_id] += hits
        for detector_id, hits in detector_hits.items():
            if hits:
                f[f'detector_hits/{detector_id}'] = np.array([hit.full_properties for hit in hits])
                f[f'detector_hits/{detector_id}'].attrs['description'] = hits[0].full_properties_description

    if visualize:
        print(f"Plotting {len(result_particles)} particles...")
        plot_particles(result_particles, plot_trajectories=track_trajectories, dimensions=3)


if __name__ == '__main__':
    def natural_number(x):
        x = int(x)
        if x < 1:
            raise argparse.ArgumentTypeError("Minimum is 1")
        return x

    parser = argparse.ArgumentParser(prog='cminject',
                                     formatter_class=argparse.MetavarTypeHelpFormatter)
    parser.add_argument('-n', help='Number of particles', type=natural_number, required=True)
    parser.add_argument('-v', help='Initial average velocity (in Z direction)', type=float, required=True)
    parser.add_argument('-o', help='Output filename for phase space (hdf5 format)', type=str, required=True)
    parser.add_argument('-f', help='Flow field filename (hdf5 format)', type=str, required=True)
    parser.add_argument('-t', help='Store trajectories?', action='store_true')
    parser.add_argument('-V', help='Visualize the particles after running?', action='store_true')
    parser.add_argument('-s', help='Run single threaded? CAUTION: Very slow, only for debugging purposes',
                        action='store_true')
    args = parser.parse_args()

    main(vz=args.v, nof_particles=args.n, flow_field=args.f, output_file=args.o,
         track_trajectories=args.t, single_threaded=args.s, visualize=args.V)


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""