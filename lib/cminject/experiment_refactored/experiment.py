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
from typing import List, Tuple, Optional, Callable
from multiprocessing import Pool, Manager
from functools import partial
import pprofile

import numpy as np
from scipy.integrate import ode

from cminject.experiment_refactored.basic import infinite_interval
from cminject.experiment_refactored.base_classes import Particle, Source, Device, Detector, ZBoundedMixin


def calculate_v_and_a(time: float, position_and_velocity: np.array,
                      particle: Particle, devices: List[Device]):
    total_acceleration = np.zeros(3, dtype=float)
    for device in devices:
        if device.is_particle_inside(particle):
            total_acceleration += device.field.calculate_acceleration(particle, time)
    return np.concatenate((position_and_velocity[3:], total_acceleration))


def is_particle_lost(particle: Particle, z_boundary: Tuple[float, float], devices: List[Device]):
    if particle.position[2] < z_boundary[0] or particle.position[2] > z_boundary[1]:
        # Particles that aren't within the experiment's total Z boundary are considered lost
        return True
    else:
        for device in devices:
            device_z_boundary = device.get_z_boundary()
            is_inside = device.is_particle_inside(particle)
            if is_inside:
                # Particles inside devices are not considered lost
                return False
            if not is_inside and device_z_boundary[0] <= particle.position[2] <= device_z_boundary[1]:
                # Particles that are outside devices *but* within their Z boundary are considered lost
                return True
        # Particles that are not inside any device but also not within any device's Z boundary are not considered
        # lost (as long as they are within the experiment's total Z boundary, see beginning of function)
        return False


def simulate_particle(particle: Particle, devices: List[Device], detectors: List[Detector],
                      t_start: float, t_end: float, dt: float,
                      z_boundary: Tuple[float, float],
                      track_trajectories: bool,
                      calc_additional: Callable[[Particle, List[Device], float], None]):
    integral = ode(calculate_v_and_a)
    integral.set_integrator('lsoda')  # TODO maybe use other integrators?
    integral.set_initial_value(np.concatenate([particle.position, particle.velocity]), t_start)
    integral.set_f_params(particle, devices)
    print(f"\tSimulating particle {particle.identifier}...")

    while integral.successful() and integral.t < t_end and not particle.lost:
        particle.time = integral.t
        # Check conditions for having lost particle.
        if not is_particle_lost(particle, z_boundary, devices):
            # If particle is not lost:
            integral_result = integral.y
            # - Store the position and velocity calculated by the integrator on the particle
            particle.position = integral_result[:3]
            particle.velocity = integral_result[3:]

            # - If present, have the calc_additional function calculate and store what it wants
            if calc_additional is not None:
                calc_additional(particle, devices, dt)

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
            print(f"\t\tParticle {particle.identifier} lost at {particle.position}. "
                  f"Tracked {len(particle.trajectory)} positions.")

    print(f"\tDone simulating particle {particle.identifier}.")
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
                 dt: float = 1.0e-5, t_start: float = 0.0, t_end: float = 1.8, track_trajectories: bool = False,
                 z_boundary: Optional[Tuple[float, float]] = None, delta_z_end: float = 0.0,
                 calc_additional: Callable[[Particle, List[Device], float], None] = None):
        # Store references
        self.devices = devices
        self.sources = sources
        self.detectors = detectors
        self.calc_additional = calc_additional

        # Store numerical data
        self.dt = dt
        self.t_start = t_start
        self.t_end = t_end
        self.track_trajectories = track_trajectories

        # Initialise list of particles from all sources
        self.particles = []
        for source in self.sources:
            source_particles = source.generate_particles(self.t_start)
            self.particles += source_particles

        # Initialise the min/max Z values amongst all detectors, and amongst all devices
        self.detectors_z_boundary = self._gather_z_boundary(self.detectors)
        self.devices_z_boundary = self._gather_z_boundary(self.devices)

        # Initialise the min/max Z boundary of the whole experiment
        if z_boundary is not None:
            # If a z_boundary was passed, use it
            self.z_boundary = z_boundary
        else:
            # Otherwise, use the global min/max from all devices and detectors
            self.z_boundary = min(self.detectors_z_boundary[0], self.devices_z_boundary[0]) - delta_z_end,\
                              max(self.detectors_z_boundary[1], self.devices_z_boundary[1]) + delta_z_end

        if self.z_boundary[0] > self.z_boundary[1]:
            exit(f"ERROR: The Z boundary of the experiment {self.z_boundary} is not sensible! The 'left' boundary "
                 f"has to be smaller than the 'right' boundary.")

    @staticmethod
    def _gather_z_boundary(z_bounded_objects: List[ZBoundedMixin]):
        max_z, min_z = infinite_interval
        for obj in z_bounded_objects:
            device_min_z, device_max_z = obj.get_z_boundary()
            min_z = min(min_z, device_min_z)
            max_z = max(max_z, device_max_z)

        return min_z, max_z

    def run_single_threaded(self, do_profiling=False):
        simulate = partial(
            (simulate_particle if not do_profiling else profiling_wrapper),
            devices=self.devices,
            t_start=self.t_start,
            t_end=self.t_end,
            dt=self.dt,
            track_trajectories=self.track_trajectories,
            detectors=self.detectors,
            z_boundary=self.z_boundary,
            calc_additional=self.calc_additional
        )
        result = map(simulate, self.particles)
        return list(result)

    def run(self, do_profiling=False, single_threaded=False):
        if single_threaded:
            return self.run_single_threaded(do_profiling=do_profiling)

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
                z_boundary=self.z_boundary,
                calc_additional=self.calc_additional
            )
            result = pool.imap_unordered(parallel_simulate, self.particles)
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