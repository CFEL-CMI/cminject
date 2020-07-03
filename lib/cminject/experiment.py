#!/usr/bin/env python3
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
import multiprocessing
import os
import random
import warnings
from typing import List, Tuple, Optional

import numpy as np
from scipy.integrate import ode

from cminject.definitions.base import ZBounded
from cminject.definitions.devices.base import Device
from cminject.definitions.detectors.base import Detector
from cminject.definitions.particles.base import Particle
from cminject.definitions.property_updaters.base import PropertyUpdater
from cminject.definitions.result_storages.base import ResultStorage
from cminject.definitions.sources.base import Source
from cminject.definitions.particle_status import *

from cminject.definitions.util import infinite_interval
from cminject.global_config import GlobalConfig, ConfigKey
from cminject.utils.args import auto_time_step
from tqdm import tqdm

# FIXME file should maybe be split into parts

PARTICLE_STATES_CONSIDERED_AS_LOST = [PARTICLE_STATUS_LOST, PARTICLE_STATUS_OUTSIDE_EXPERIMENT]

DEVICES, DETECTORS, PROPERTY_UPDATERS, Z_BOUNDARY, NUMBER_OF_DIMENSIONS, BASE_SEED = [None] * 6
TIME_START, TIME_END, TIME_STEP = [None] * 3
Z_LOWER, Z_UPPER = None, None


class ImpossiblePropagationException(Exception):
    pass


def get_particle_simulation_state(particle_position: np.array, time: float,
                                  z_boundary: Tuple[float, float], devices: List[Device]) -> int:
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
        return PARTICLE_STATUS_OUTSIDE_EXPERIMENT
    else:
        for device in devices:
            device_z_boundary = device.z_boundary
            is_inside = device.is_particle_inside(particle_position, time)
            if is_inside:
                # Particles inside devices are not considered lost
                return PARTICLE_STATUS_OK
            if not is_inside and device_z_boundary[0] <= particle_z_pos <= device_z_boundary[1]:
                # Particles that are outside devices *but* within their Z boundary are considered lost
                return PARTICLE_STATUS_LOST
        # Particles that are not inside any device but also not within any device's Z boundary are not considered
        # lost (as long as they are within the experiment's total Z boundary, see beginning of function)
        return PARTICLE_STATUS_BETWEEN_DEVICES


def postprocess_integration_step(particle: Particle, integral: ode, position_changed: bool) -> None:
    """
    Does post-processing after an integration step, which consists of:
      - running all property updaters, possibly re-instantiating the integrator when they change the particle position
      - running all detectors, asking them to try and detect the particle

    :param particle: The particle which experienced an integration step
    :param integral: The integral (scipy.integrate.ode) instance.
    :param position_changed:
        Whether the position of the particle had been changed relative to the integral result,
        _before_ this function was called.
    """
    # Run all property updaters for the particle with the current integral time
    for property_updater in PROPERTY_UPDATERS:
        position_changed |= property_updater.update(particle, integral.t)
    # If some property updater changed the position, reset the integrator's initial value
    # NOTE: doing this is slow and should be reserved for cases where motion cannot be modeled differently
    if position_changed:
        integral.set_initial_value(particle.position, integral.t)

    # Have each detector try to detect the particle
    for detector in DETECTORS:
        detector.try_to_detect(particle)


def postprocess_particle(particle: Particle, integral: ode, t_end: float) -> None:
    # Run the updaters and detectors one last time after simulation completed
    for property_updater in PROPERTY_UPDATERS:
        property_updater.update(particle, integral.t)
    for detector in DETECTORS:
        detector.try_to_detect(particle)

    # Log some additional info if loglevel is low enough
    if logging.root.level <= logging.INFO:
        if particle.lost:
            z_position = particle.spatial_position[NUMBER_OF_DIMENSIONS - 1]
            if Z_BOUNDARY[0] <= z_position <= Z_BOUNDARY[1]:
                reason = 'hit boundary within the experiment'
            else:
                reason = f'left experiment Z boundary (at {z_position:.2g})'
        elif integral.t >= t_end:
            reason = 'whole timespan simulated'
        else:
            reason = 'unknown reason'

        logging.debug(f"Last particle position: {particle.position}")
        logging.info(f"\tDone simulating particle {particle.identifier}: {reason}.")


def propagate_field_free(particle: Particle) -> None:
    """
    Propagates a particle assuming it is in a field-free region, up until the next valid point in Z direction.
      - "Z direction" is always the last dimension, assumed to be the main propagation axis.
      - "the next valid point" is calculated from the z_boundary properties of all devices and detectors.
      - Propagation is done according to a simple rule-of-three calculation.

    :param particle: The Particle instance to propagate.
    :raises ImpossiblePropagationException: if asked to propagate the particle in an impossible way, e.g.
        propagating until the next point in positive Z when the particle's v_Z is negative.
    """
    z, vz = particle.spatial_position[NUMBER_OF_DIMENSIONS - 1], particle.velocity[NUMBER_OF_DIMENSIONS - 1]
    if vz > 0:
        if z <= Z_LOWER[0]:
            raise ImpossiblePropagationException(
                f"Cannot propagate field-free from {z} until {Z_LOWER[0]} with v_z > 0."
            )
        idx = np.searchsorted(Z_LOWER, z, side='right')
        next_z = Z_LOWER[idx]
    elif vz < 0:
        if z >= Z_UPPER[-1]:
            raise ImpossiblePropagationException(
                f"Cannot propagate field-free from {z} until {Z_UPPER[-1]} with v_z < 0."
            )
        idx = len(Z_UPPER)-1 - np.searchsorted(-Z_UPPER[::-1], -z, side='right')  # TODO hmm... maybe preprocess
        next_z = Z_UPPER[idx]
    else:  # vz == 0
        raise ImpossiblePropagationException("Cannot propagate field-free with v_z = 0.")

    prevpos = np.copy(particle.position)
    delta_t = abs(next_z - z) / abs(vz)
    particle.position = np.concatenate([
        particle.spatial_position + delta_t * particle.velocity,
        particle.velocity
    ])
    logging.info(f"Field-free propagation from {z} until {next_z} with velocities={particle.velocity} and "
                 f"delta_t={delta_t}.\nBefore: {prevpos}. After: {particle.position}.")


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
    np.random.seed(BASE_SEED + int(particle.identifier))  # reseed RandomState from the base seed and particle ID

    # Construct integrals
    integral = ode(spatial_derivatives)
    integral.set_integrator('lsoda', nsteps=3000)
    integral.set_initial_value(particle.position, TIME_START)
    integral.set_f_params(particle, DEVICES, NUMBER_OF_DIMENSIONS)
    logging.info(f"\tSimulating particle {particle.identifier}...")

    # TODO:
    # when loop was broken due to the integration not being successful, we need to somehow continue to simulate
    # the particle another way, e.g. with another integrator (how? which?)
    while integral.successful() and integral.t < TIME_END and not particle.lost:
        # Store the position and velocity calculated by the integrator on the particle
        integral_result = integral.y
        particle.position = integral_result
        particle.time_of_flight = integral.t
        postprocess_integration_step(particle, integral, position_changed=False)

        # Check what 'simulation state' the particle is in wrt. the devices, experiment Z boundary, and integral time
        simulation_state = get_particle_simulation_state(
            particle.position[:NUMBER_OF_DIMENSIONS], integral.t, Z_BOUNDARY, DEVICES
        )

        if simulation_state in PARTICLE_STATES_CONSIDERED_AS_LOST:
            # Store that particle is lost and exit the loop
            particle.lost = True
            break
        elif simulation_state is PARTICLE_STATUS_BETWEEN_DEVICES:
            try:
                propagate_field_free(particle)
            except ImpossiblePropagationException as e:
                logging.error(f"Impossible propagation for particle {particle.identifier}: {e}")
                particle.lost = True
                break
            # Since we just "simulated" an integration step, we have to run the postprocessing again
            postprocess_integration_step(particle, integral, position_changed=True)

        # Propagate by integrating until the next time step
        integral.integrate(integral.t + TIME_STEP)

    postprocess_particle(particle, integral, TIME_END)
    return particle


class Experiment:
    def __init__(self, number_of_dimensions: int, time_interval: Tuple[float, float], time_step: float = 1e-5,
                 z_boundary: Optional[Tuple[float, float]] = None, random_seed=None,
                 result_storage: Optional[ResultStorage] = None):
        """
        Construct an Experiment to run.

        :param devices: The list of Devices in the experimental setup.
        :param sources: The list of Sources that generate particles. Note that for both simulation and storage reasons,
            having different sources generate different types of particles for one Experiment is not allowed.
        :param detectors: The list of Detectors in the experimental setup.
        :param property_updaters: The list of PropertyUpdaters used to update different properties on the particles
            in each step.
        :param result_storage: The ResultStorage to use for the experiment's results. Can be None, in which case
            the results will not be stored and only returned.
        :param time_interval: The time interval to run the experiment in, as a 2-tuple of (`t_start`, `t_end`).
        :param time_step: The time step to use. If None is passed, a time step is picked based on the initial particle
            velocities.
        :param z_boundary: (Optional) If passed, the Z boundary of the entire experimental setup. Particles will not
            be simulated outside of this boundary. If this argument is not passed, the minimal interval enclosing all
            Devices and Detectors along the Z axis will be calculated and used as the Z boundary.
        """
        self.devices = []
        self.sources = []
        self.detectors = []
        self.property_updaters = []
        self.result_storage = result_storage
        self.random_seed = random_seed
        self.time_interval = time_interval
        self.time_step = time_step
        self.z_boundary = z_boundary
        self.loglevel = 'warning'

        # TODO do we need to subscribe?
        self.number_of_dimensions = number_of_dimensions
        GlobalConfig().set(ConfigKey.NUMBER_OF_DIMENSIONS, number_of_dimensions)
        GlobalConfig().set(ConfigKey.TIME_STEP, time_step)

    def add_device(self, device: Device):
        """
        Adds a Device to this experiment.

        :param device: The Device instance to add.
        """
        self.devices.append(device)
        return self

    def add_source(self, source: Source):
        """
        Adds a Source to this experiment.

        :param source: The Source instance to add.
        """
        self.sources.append(source)
        return self

    def add_detector(self, detector: Detector):
        """
        Adds a Detector to this experiment.

        :param detector: The Detector instance to add.
        """
        self.detectors.append(detector)
        return self

    def add_property_updater(self, property_updater: PropertyUpdater):
        """
        Adds a PropertyUpdater to this experiment.

        :param property_updater: The PropertyUpdater instance to add.
        """
        self.property_updaters.append(property_updater)
        return self

    @property
    def time_interval(self):
        return self.__time_interval

    @time_interval.setter
    def time_interval(self, val):
        if not isinstance(val, tuple) or not len(val) == 2:
            raise ValueError("time_interval must be a tuple of two floats (t_begin, t_end)!")
        time_begin, time_end = val
        if time_begin >= time_end:
            errmsg = f"time_begin must be less than time_end, but time_begin={time_begin}, time_end={time_end}!"
            raise ValueError(errmsg)
        self.__time_interval = (time_begin, time_end)

    @property
    def time_step(self):
        return self.__time_step

    @time_step.setter
    def time_step(self, val):
        if val <= 0:
            raise ValueError(f"time_step must be > 0, but {val} was given!")
        self.__time_step = val

    @property
    def z_boundary(self):
        return self.__z_boundary

    @z_boundary.setter
    def z_boundary(self, val):
        if val is None:
            self.__z_boundary = None
            return

        if not isinstance(val, tuple) or not len(val) == 2:
            raise ValueError("z_boundary must be a tuple of two floats (z_start, z_end)!")
        z_start, z_end = val
        if z_start >= z_end:
            errmsg = f"z_start must be less than z_end, but z_start={z_start}, z_end={z_end}!"
            raise ValueError(errmsg)
        self.__z_boundary = (z_start, z_end)

    def run(self,
            single_threaded: bool = False, chunksize: int = None, processes: int = os.cpu_count(),
            loglevel: str = "warning",
            progressbar: bool = False) -> List[Particle]:
        """
        Run the Experiment. A list of resulting Particle instances is returned.

        :param single_threaded: Whether to run this experiment on a single thread, i.e. without parallelization. False
            by default: It's mostly only useful to pass True for developers, as many debuggers, profilers, etc. are not
            able to deal well with different processes/threads. To run the experiment and get results quickly, leave
            this on False.
        :param chunksize: The chunk size for pool.imap_unordered. None by default, which equates to 1. This parameter
            has no effect if single_threaded is True.
        :param processes: The number of processes to use for multiprocessing. Handed directly to multiprocessing.Pool(),
            so check the documentation there for any further info. Equal to the number of CPU cores by default.
        :param loglevel: The loglevel to run the program with. One of {DEBUG, INFO, WARNING, ERROR, CRITICAL} --
            see Python's builtin logging module for more explanation regarding these levels.
        :param progressbar: Show simulation progress in a readable progress bar (by tqdm).
        :return: A list of resulting Particle instances. Things like detector hits and trajectories should be stored on
            them and can be read off each Particle.
        """
        self.loglevel = loglevel

        self._prepare()
        if single_threaded:
            logging.info("Running single-threaded.")
            iterator = map(simulate_particle, self.particles)
            if progressbar:
                iterator = tqdm(iterator, total=len(self.particles))
            particles = list(iterator)
        else:
            logging.info(f"Running in parallel using {processes} processes with chunksize {chunksize}.")
            pool = multiprocessing.Pool(processes=processes, initializer=self._prepare_thread, initargs=())
            try:
                iterator = pool.imap_unordered(simulate_particle, self.particles, chunksize=chunksize)
                if progressbar:
                    iterator = tqdm(iterator, total=len(self.particles))
                particles = list(iterator)
            finally:
                pool.close()
                pool.join()

        if self.result_storage is not None:
            with self.result_storage as s:
                s.store_results(particles)

        return particles

    def _prepare(self):
        # Set values on global config, thereby cascading them throughout the setup to all subscribed objects
        conf = GlobalConfig()
        conf.set(ConfigKey.NUMBER_OF_DIMENSIONS, self.number_of_dimensions)
        conf.set(ConfigKey.TIME_STEP, self.time_step)

        # Initialise list of particles from all sources
        particle_types = set()  # Remember a set of particle types to raise an error if there are multiple types
        self.particles = []
        for source in self.sources:
            source_particles = source.generate_particles(self.time_interval[0])  # initialize with t_0
            if source_particles:
                particle_types.add(type(source_particles[0]))
            self.particles += source_particles
        if len(particle_types) > 1:
            raise TypeError("Cannot simulate particles of multiple types at the same time!")

        # Initialize random state
        if self.random_seed is None:
            seed = np.random.randint(2 ** 31)
            self.random_seed = seed
        np.random.seed(self.random_seed)

        # Initialise the min/max Z values amongst all detectors and all devices
        self.detectors_z_boundary = self._gather_z_boundary(self.detectors)
        self.devices_z_boundary = self._gather_z_boundary(self.devices)
        self.z_lower, self.z_upper = self._gather_z_lower_upper([
            x.z_boundary for x in self.detectors + self.devices
        ])

        # Initialise the min/max Z boundary of the whole experiment
        if self.z_boundary is None:
            # Use the global min/max from all devices and detectors, -/+ the end Z delta
            self.z_boundary = min(self.detectors_z_boundary[0], self.devices_z_boundary[0]), \
                              max(self.detectors_z_boundary[1], self.devices_z_boundary[1])

        self._prepare_thread()

    def _prepare_thread(self):
        self._set_globals()

        current_process = multiprocessing.current_process()
        if isinstance(current_process, (multiprocessing.context.ForkProcess, multiprocessing.context.SpawnProcess)):
            mp_id = current_process.name
            log_format = '%(levelname)s:[%(filename)s/%(funcName)s@' + mp_id + '] %(message)s'
        else:
            log_format = '%(levelname)s:[%(filename)s/%(funcName)s] %(message)s'

        # we do the logging initialization in this function so it is done in each runner thread as well
        logging.basicConfig(format=log_format, level=getattr(logging, self.loglevel.upper()))

    def _set_globals(self):
        global DEVICES, DETECTORS, PROPERTY_UPDATERS, NUMBER_OF_DIMENSIONS, BASE_SEED
        global TIME_START, TIME_END, TIME_STEP
        global Z_BOUNDARY, Z_LOWER, Z_UPPER
        DEVICES = self.devices
        DETECTORS = self.detectors
        PROPERTY_UPDATERS = self.property_updaters
        TIME_START, TIME_END = self.time_interval
        TIME_STEP = self.time_step
        Z_BOUNDARY = self.z_boundary
        Z_LOWER, Z_UPPER = self.z_lower, self.z_upper
        NUMBER_OF_DIMENSIONS = self.number_of_dimensions
        BASE_SEED = self.random_seed

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

    @staticmethod
    def _gather_z_lower_upper(z_boundaries):
        """
        An internal method to gather the lower and upper bounds of all intervals, grouped into an array
        of lower bounds and another array of upper bounds.

        :param z_bounded_objects: A list of Z bounded objects.
        :return: A 2-tuple of (min_z, max_z).
        """
        arr = np.array(z_boundaries)
        lower = np.array(list(sorted(arr[:, 0])))
        upper = np.array(list(sorted(arr[:, 1], reverse=True)))
        logging.debug(f"Lower Z boundaries: {lower}")
        logging.debug(f"Upper Z boundaries: {upper}")
        return lower, upper

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
