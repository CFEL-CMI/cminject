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
import logging
from functools import partial
from multiprocessing import Pool
from typing import List, Tuple, Optional

import numpy as np
from cminject.definitions.base import Particle, Source, Device, Detector, ZBounded, \
    PropertyUpdater, infinite_interval
from scipy.integrate import ode


def spatial_derivatives(time: float, position_and_velocity: np.array,
                        particle: Particle, devices: List[Device], number_of_dimensions: int) -> np.array:
    """
    Calculates the derivatives of position and velocity, which will then be integrated over.
    "n" as used below is the number of the simulation's spatial dimensions.

    :param time: The current simulation time in seconds
    :param position_and_velocity: A (2n,) numpy array containing spatial position and velocity of the particle.
    :param particle: The Particle instance.
    :param devices: The list of devices that (might) affect this particle.
    :param number_of_dimensions: The number of dimensions of the space the particle moves in.
    :return: A (2n,) numpy array containing velocity and acceleration of the particle.
    """
    particle.position = position_and_velocity
    total_acceleration = np.zeros(number_of_dimensions, dtype=float)
    for device in devices:
        if device.is_particle_inside(position_and_velocity[:number_of_dimensions], time):
            total_acceleration += device.calculate_acceleration(particle, time)
    return np.concatenate((position_and_velocity[number_of_dimensions:], total_acceleration))


def is_particle_lost(particle_position: np.array, time: float,
                     z_boundary: Tuple[float, float], devices: List[Device]):
    """
    Decides whether a particle should be considered lost.

    :param particle_position: The particle's position.
    :param time: The current simulation time in seconds
    :param z_boundary: The total Z boundary to consider
    :param devices: The list of devices that the particle might potentially be in
    :return:
    """
    particle_z_pos = particle_position[-1]
    if particle_z_pos < z_boundary[0] or particle_z_pos > z_boundary[1]:
        # Particles that aren't within the experiment's total Z boundary are considered lost
        return True
    else:
        for device in devices:
            device_z_boundary = device.z_boundary
            is_inside = device.is_particle_inside(particle_position, time)
            if is_inside:
                # Particles inside devices are not considered lost
                return False
            if not is_inside and device_z_boundary[0] <= particle_z_pos <= device_z_boundary[1]:
                # Particles that are outside devices *but* within their Z boundary are considered lost
                return True
        # Particles that are not inside any device but also not within any device's Z boundary are not considered
        # lost (as long as they are within the experiment's total Z boundary, see beginning of function)
        return False


DEVICES, DETECTORS, PROPERTY_UPDATERS, TIME_INTERVAL, Z_BOUNDARY, NUMBER_OF_DIMENSIONS, BASE_SEED = [None] * 7


def simulate_particle(particle: Particle) -> Particle:
    """
    Simulates the flight path of a single particle, with a list of devices that can affect the particle,
    a list of detectors that can detect the particle, a list of property updaters that can store additional results
    on the particle, a time interval to simulate the particle for, and a Z boundary to simulate the particle within.

    :param particle: The particle instance.
    :return: A modified version of the particle instance after the simulation has ended.

    .. note::
        There are implicit parameters, as global variables. This is done to avoid these (unchanging) objects being
        serialized and deserialized for each particle that is being simulated, see this blog post:
        https://confluence.desy.de/display/PARTI/2019/08/08/CMInject%3A+Performance+considerations

        For more detailed info on this issue in general, see:
        https://thelaziestprogrammer.com/python/a-multiprocessing-pool-pickle

        For the way this is worked around using global variables and a Pool initializer, see:
        https://stackoverflow.com/a/10118250/3090225

        This greatly improves performance, as we're not doing unnecessary I/O during which we can't do any actual
        calculations. Below is the list of global variables that used to be parameters bound via functools.partial,
        but are now global variables that are used within the function:

        - DEVICES: The list of devices.
        - DETECTORS: The list of detectors.
        - PROPERTY_UPDATERS: The list of property_updaters.
        - TIME_INTERVAL: The time interval of the simulation.
        - Z_BOUNDARY: The Z boundary of the simulation.
        - NUMBER_OF_DIMENSIONS: The number of dimensions of the space the particle moves in.
        - BASE_SEED: The "base seed", i.e. the random seed the experiment was defined with, to derive a local
          random seed from.
    """
    t_start, t_end, dt = TIME_INTERVAL
    np.random.seed(BASE_SEED + int(particle.identifier))  # reseed RandomState from the base seed and particle ID

    # Run all property updaters once with the start time
    for property_updater in PROPERTY_UPDATERS:
        property_updater.update(particle, t_start)

    # Run all detectors once
    for detector in DETECTORS:
        detector.try_to_detect(particle)

    # Construct integrals
    integral = ode(spatial_derivatives)
    integral.set_integrator('lsoda')  # TODO maybe use other integrators?
    integral.set_initial_value(particle.position, t_start)
    integral.set_f_params(particle, DEVICES, NUMBER_OF_DIMENSIONS)
    logging.info(f"\tSimulating particle {particle.identifier}...")

    # TODO:
    """
    when loop was broken due to the integration not being successful, we need to somehow continue to simulate
    the particle another way, e.g. with another integrator (how? which?)
    """
    while integral.successful() and integral.t < t_end and not particle.lost:
        # Check conditions for having lost particle.
        if not is_particle_lost(particle.position[:NUMBER_OF_DIMENSIONS], integral.t, Z_BOUNDARY, DEVICES):
            # If particle is not lost:
            # - Store the position and velocity calculated by the integrator on the particle
            integral_result = integral.y
            particle.position = integral_result
            particle.time_of_flight = integral.t

            # - Run all property updaters for the particle with the current integral time
            changed_position = False
            for property_updater in PROPERTY_UPDATERS:
                changed_position |= property_updater.update(particle, integral.t)

            # - Have each detector try to detect the particle
            for detector in DETECTORS:
                detector.try_to_detect(particle)

            # Propagate by integrating until the next time step
            if changed_position:
                integral.set_initial_value(particle.position, integral.t)
            integral.integrate(integral.t + dt)
        else:
            # If particle is lost, store this and (implicitly) break the loop
            particle.lost = True

    if particle.lost:
        if Z_BOUNDARY[0] <= particle.spatial_position[NUMBER_OF_DIMENSIONS - 1] <= Z_BOUNDARY[1]:
            reason = 'hit boundary within the experiment'
        else:
            reason = 'left experiment Z boundary'
    elif integral.t >= t_end:
        reason = 'whole timespan simulated'
    else:
        reason = 'unknown reason'

    logging.info(f"\tDone simulating particle {particle.identifier}: {reason}.")
    return particle


class Experiment:
    """"""

    def __init__(self, devices: List[Device], sources: List[Source], detectors: List[Detector],
                 property_updaters: List[PropertyUpdater] = None,
                 time_interval: Tuple[float, float, float] = (0.0, 1.8, 1.0e-6),
                 z_boundary: Optional[Tuple[float, float]] = None, delta_z_end: float = 0.0,
                 number_of_dimensions: int = 3, seed=None):
        """
        Construct an Experiment to run.

        :param devices: The list of Devices in the experimental setup.
        :param sources: The list of Sources that generate particles. Note that for both simulation and storage reasons,
            having different sources generate different types of particles for one Experiment is not allowed.
        :param detectors: The list of Detectors in the experimental setup.
        :param property_updaters: The list of PropertyUpdaters used to update different properties on the particles
            in each step.
        :param time_interval: The time interval to run the experiment in, in the shape of (`t_start`, `t_end`, `dt`).
            The simulation will run in the time interval [`t_start`, `t_end`] with time step `dt`.
        :param z_boundary: (Optional) If passed, the Z boundary of the entire experimental setup. Particles will not
            be simulated outside of this boundary. If not passed, the minimal interval enclosing all Devices and
            Detectors along the Z axis will be calculated and used as the Z boundary.
        :param delta_z_end: An additional amount to add on both ends of the Z boundary. Mostly sensible to use when
            not passing an explicit `z_boundary`, and when you expect the automatically generated minimal boundary to be
            too small, e.g. if some Source will generate particles that lie slightly outside the automatically generated
            Z boundary.
        """
        # Store references
        self.devices = devices
        self.sources = sources
        self.detectors = detectors
        self.property_updaters = property_updaters

        # Store the number of dimensions and set it on all relevant objects
        self.number_of_dimensions = number_of_dimensions
        for dimensional_object in (self.devices + self.sources + self.detectors + self.property_updaters):
            dimensional_object.set_number_of_dimensions(self.number_of_dimensions)

        # Store the time interval and the seed, use the seed to seed numpy's RandomState
        self.time_interval = time_interval
        if seed is None:
            seed = np.random.randint(2 ** 31)
        self.seed = seed
        np.random.seed(self.seed)

        # Initialise list of particles from all sources
        particle_types = set()  # Remember a set of particle types to raise an error if there are multiple
        self.particles = []
        for source in self.sources:
            source_particles = source.generate_particles(self.time_interval[0])  # 0 is t_start
            if source_particles:
                # FIXME assuming a source only generates one particle type
                particle_types.add(type(source_particles[0]))

            # Set the number of dimensions for each particle.
            for source_particle in source_particles:
                source_particle.set_number_of_dimensions(number_of_dimensions)

            self.particles += source_particles
        if len(particle_types) > 1:
            raise TypeError("Cannot simulate particles of multiple types at the same time!")

        # Initialise the min/max Z values amongst all detectors, and amongst all devices
        self.detectors_z_boundary = self._gather_z_boundary(self.detectors)
        self.devices_z_boundary = self._gather_z_boundary(self.devices)

        # Initialise the min/max Z boundary of the whole experiment
        if z_boundary is not None:
            # If a z_boundary was passed, use it
            self.z_boundary = z_boundary
        else:
            # Otherwise, use the global min/max from all devices and detectors, -/+ the end Z delta
            self.z_boundary = min(self.detectors_z_boundary[0], self.devices_z_boundary[0]) - delta_z_end, \
                              max(self.detectors_z_boundary[1], self.devices_z_boundary[1]) + delta_z_end

        if self.z_boundary[0] > self.z_boundary[1]:
            exit(f"ERROR: The Z boundary {self.z_boundary} of the experiment is not sensible! The 'left' boundary "
                 f"has to be smaller than the 'right' boundary.")

    def run(self, single_threaded=False) -> List[Particle]:
        """
        Run the Experiment. A list of resulting Particle instances is returned.

        :param single_threaded: Whether to run this experiment on a single thread, i.e. without parallelization. False
            by default: It's mostly only useful to pass True for developers, as many debuggers, profilers, etc. are not
            able to deal well with different processes/threads. To run the experiment and get results quickly, leave
            this on False.
        :return: A list of resulting Particle instances. Things like detector hits and trajectories should be stored on
            them and can be read off each Particle.
        """

        if single_threaded:
            self._initialize_globals()
            particles = list(map(simulate_particle, self.particles))
        else:
            pool = Pool(initializer=self._initialize_globals, initargs=())
            try:
                particles = list(pool.imap_unordered(simulate_particle, self.particles))
            finally:
                pool.close()
                pool.join()

        return particles

    @staticmethod
    def _gather_z_boundary(z_bounded_objects: List[ZBounded]) -> Tuple[float, float]:
        """
        An internal method to gather the minimal enclosing Z boundary for a list of Z bounded objects.

        :param z_bounded_objects: A list of Z bounded objects.
        :return: A 2-tuple of (min_z, max_z).
        """
        max_z, min_z = infinite_interval
        for obj in z_bounded_objects:
            obj_min_z, obj_max_z = obj.z_boundary
            min_z = min(min_z, obj_min_z)
            max_z = max(max_z, obj_max_z)
        return min_z, max_z

    def _initialize_globals(self):
        global DEVICES, DETECTORS, PROPERTY_UPDATERS, TIME_INTERVAL, Z_BOUNDARY, NUMBER_OF_DIMENSIONS, BASE_SEED
        DEVICES = self.devices
        DETECTORS = self.detectors
        PROPERTY_UPDATERS = self.property_updaters
        TIME_INTERVAL = self.time_interval
        Z_BOUNDARY = self.z_boundary
        NUMBER_OF_DIMENSIONS = self.number_of_dimensions
        BASE_SEED = self.seed


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
