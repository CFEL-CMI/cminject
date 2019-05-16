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

from typing import List, Tuple
from multiprocessing import Pool, Manager
from functools import partial
import pprofile

import numpy as np
from scipy.integrate import ode

from .base_classes import Particle, Source, Device, Detector


def calculate_v_and_a(time: float, position_and_velocity: np.array,
                      particle: Particle, devices: List[Device]):
    total_acceleration = np.zeros(3, dtype=float)
    for device in devices:
        if not device.is_particle_inside(particle):
            continue
        total_acceleration += device.field.calculate_acceleration(particle, time)
    return np.concatenate((position_and_velocity[3:], total_acceleration))


def is_particle_lost(particle: Particle, z_boundary: Tuple[float, float], devices: List[Device]):
    if z_boundary[0] <= particle.position[2] <= z_boundary[1]:
        return False
    else:
        return not any(device.boundary.is_particle_inside(particle) for device in devices)


def simulate_particle(particle: Particle, devices: List[Device], detectors: List[Detector],
                      t_start: float, t_end: float, dt: float,
                      z_boundary: Tuple[float, float],
                      track_trajectories: bool):
    integral = ode(calculate_v_and_a)
    integral.set_integrator('lsoda')  # TODO maybe use other integrators?
    integral.set_initial_value(np.concatenate([particle.position, particle.velocity]), t_start)
    integral.set_f_params(particle, devices)

    while integral.successful() and integral.t < t_end and not particle.lost:
        # Check conditions for having lost particle.
        if not is_particle_lost(particle, z_boundary, devices):
            # If particle is not lost:
            integral_result = integral.y
            # - Store the position and velocity calculated by the integrator on the particle
            particle.position = integral_result[:3]
            particle.velocity = integral_result[3:]

            # - Have each detector store a hit position if there was a hit
            for detector in detectors:
                if detector.has_reached_detector(particle):
                    hit_position = detector.get_hit_position(particle)
                    particle.store_hit(detector.identifier, hit_position)

            # - If we want to track trajectories, store them
            if track_trajectories:
                particle.trajectory.append(
                    np.concatenate([[integral.t], integral_result[0:6]])
                )

            integral.integrate(integral.t + dt)
        else:
            # If particle is lost, store this and (implicitly) break the loop
            particle.lost = True

    return particle


def profiling_wrapper(particle: Particle, devices: List[Device], detectors: List[Detector],
                      t_start: float, t_end: float, dt: float,
                      z_boundary: Tuple[float, float],
                      track_trajectories: bool):
    result = None
    prof = pprofile.Profile()
    with prof():
        result = simulate_particle(particle, devices, detectors,
                                   t_start, t_end, dt,
                                   z_boundary, track_trajectories)
    prof.dump_stats("cachegrind.out.%d" % particle.identifier)
    return result


class Experiment:
    def __init__(self, devices: List[Device], sources: List[Source], detectors: List[Detector],
                 dt: float = 1.0e-5, t_start=0.0, t_end=1.8, track_trajectories=False):
        # Store references
        self.devices = devices
        self.sources = sources
        self.detectors = detectors

        # Store numerical data
        self.dt = dt
        self.t_start = t_start
        self.t_end = t_end
        self.track_trajectories = track_trajectories

        # Initialise list of particles from all sources
        self.particles = []
        for source in self.sources:
            source_particles = source.generate_particles()
            self.particles += source_particles

        # Initialise the min/max Z boundary of the whole experiment
        min_z = float('inf')
        max_z = float('-inf')
        for device in self.devices:
            device_min_z, device_max_z = device.get_z_boundary()
            min_z = min(min_z, device_min_z)
            max_z = max(max_z, device_max_z)
            self.z_boundary = (min_z, max_z)

    def run(self, do_profiling=False):
        print("Running experiment...")
        pool = Pool()
        try:
            parallel_simulate = partial(
                (simulate_particle if not do_profiling else profiling_wrapper),
                devices=self.devices,
                t_start=self.t_start,
                t_end=self.t_end,
                dt=self.dt,
                track_trajectories=self.track_trajectories,
                detectors=self.detectors,
                z_boundary=self.z_boundary
            )
            result = map(parallel_simulate, self.particles)
        except Exception as e:
            raise e
        finally:
            if do_profiling:
                prof = pprofile.Profile()
                with prof():
                    pool.close()
                    pool.join()
                prof.dump_stats("cachegrind.out.main")
            else:
                pool.close()
                pool.join()

        return list(result)


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""