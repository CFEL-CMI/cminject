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

"""
The most fundamental of definitions, for deriving other base classes from.
"""

__all__ = ['ZBounded', 'Particle', 'Boundary', 'Detector', 'Device', 'Field', 'Action', 'ResultStorage',
           'Setup', 'Source', 'ParticleDetectorHit']

import argparse
import logging
from abc import ABC, abstractmethod
from typing import Tuple, Optional, List, Any, Dict, Iterable, Union

import numpy as np

from cminject.utils.args import SetupArgumentParser
from cminject.utils.global_config import ConfigSubscriber, GlobalConfig, ConfigKey
from cminject.utils.perf import cached_property


class ZBounded(ABC):
    """
    A mixin for objects in an experiment setup that are bounded in the Z direction.
    Classes deriving from this mixin will have to implement the get_z_boundary method.

    Objects implementing this mixin should be used to determine the minimum/maximum extent of the entire experiment
    setup in the Z direction, and then allows an Experiment to make a fast preliminary calculation about whether it
    can consider a particle lost based on this.
    """
    @property
    @abstractmethod
    def z_boundary(self) -> Tuple[float, float]:
        """
        Returns the Z boundary of this Z-bounded object.

        :return: A tuple of floats, the first entry being z_min, the second being z_max.
        """
        pass


class Particle(ABC):
    """
    Describes a particle whose trajectory we want to simulate.
    It is first and foremost a data container, and it and its subclasses should be written and used as such.

    It can be read by any part of the code, but should only be written to by instances of
    :class:`cminject.experiment.Experiment` and :class:`cminject.base.Action` (or instances of their subclasses).

    This class declares a few basic properties that we expect every particle to have:

      - identifier (a unique integer id),
      - lost (a flag storing whether the particle is considered lost)
      - position (the particle's position)
      - velocity (the particle's velocity)
      - mass (the particle's mass)
      - trajectory (a list describing points in the particle's path)
      - detector_hits (a dict mapping detector identifiers to hit lists)
    """
    def __init__(self, identifier: int, start_time: float, position: np.array, velocity: np.array, *args, **kwargs):
        """
        The constructor for Particle.

        :param identifier: The unique identifier the particle should have.
        :param start_time: The time the particle starts being simulated at.
        :param position: The phase space position (n-D position and velocity) the particle starts at.
            A (2*n,) dimensional vector, where n is the number of spatial dimensions.
        """
        self.identifier: int = identifier
        self.lost: bool = False
        self.time: float = start_time

        assert len(position) == len(velocity),\
            "Position and velocity must have same length but have shapes %s and %s" % (position.shape, velocity.shape)
        self.number_of_dimensions = len(position)
        self.phase_space_position: np.array = np.concatenate((position, velocity))  # this copies implicitly
        self.trajectory: List[np.array] = []
        self.detector_hits: Dict[str, List['ParticleDetectorHit']] = {}
        self._initial_tracked_properties: np.array = self.as_array('tracked')

    @property
    def position(self):
        """The spatial position of the particle."""
        return self.phase_space_position[:self.number_of_dimensions]

    @position.setter
    def position(self, value):
        self.phase_space_position[:self.number_of_dimensions] = value

    @property
    def velocity(self):
        """The spatial velocity of the particle."""
        return self.phase_space_position[self.number_of_dimensions:]

    @velocity.setter
    def velocity(self, value):
        self.phase_space_position[self.number_of_dimensions:] = value

    @property
    def initial_tracked_properties(self):
        """An array representing the state of all tracked properties of this particle, at its time of creation."""
        return self._initial_tracked_properties

    @cached_property
    @abstractmethod
    def mass(self) -> float:
        """
        The mass of the particle.

        .. warning::
          Please ensure that you also use the @cached_property decorator when overriding this property, so that the mass
          calculation is done exactly once and the result is stored, for (time) efficiency.
        """
        pass

    @cached_property
    def tracked_properties(self):
        """
        The definition of the particle's tracked properties, i.e. properties that are changing along the trajectory
        of the particle. Must return a list in `NumPy structured dtype format
        <https://numpy.org/doc/stable/user/basics.rec.html#structured-datatype-creation>`_).
        """
        return [
            ('position', (np.float64, self.number_of_dimensions)),
            ('velocity', (np.float64, self.number_of_dimensions)),
            ('time', np.float64)
        ]

    @cached_property
    def constant_properties(self):
        """
        The definition of the particle's constant properties, i.e. properties that do *NOT* change along the
        trajectory of the particle. Must return a list in `NumPy structured dtype format
        <https://numpy.org/doc/stable/user/basics.rec.html#structured-datatype-creation>`_).
        """
        return [
            ('identifier', np.int32),
            ('mass', np.float64)
        ]

    def _get_props(self, which):
        if which == 'all':
            return self.tracked_properties + self.constant_properties
        elif which == 'tracked':
            return self.tracked_properties
        elif which == 'constant':
            return self.constant_properties

    def as_dtype(self, which) -> np.array:
        """
        Returns the dtype that the array returned by as_array() will have, when called with the same parameters.

        :param which: The part of this particle's state to return. Refer to the docstring about as_array.
        """
        return np.dtype(self._get_props(which))

    def as_array(self, which) -> np.array:
        """
        Returns a NumPy array representation of the current state of this particle. Either the tracked properties,
        or the constant properties, or all (tracked + constant) can be returned as a single NumPy array with a
        structured NumPy datatype.

        :param which: The part of this particle's state to return as an array. Must be one of
         ['tracked', 'constant', 'all']:

           - 'tracked' will refer to the implementation of tracked_properties,
           - 'constant' will refer to the implementation of constant_properties,
           - 'all' will refer to both and concatenate them in order (tracked + constant).
        """
        assert which in ['tracked', 'constant', 'all'],\
            f"which was '{which}' but must be 'tracked', 'constant', or 'all'"
        props = self._get_props(which)
        return np.array([
            tuple(getattr(self, prop[0]) for prop in props)
        ], dtype=self.as_dtype(which))

    @property
    def reached_any_detector(self) -> bool:
        """
        Whether this particle has ever reached any detector.
        """
        return not not self.detector_hits  # `not not` to convert to boolean, to not leak data here

    def __str__(self):
        return f"<{self.__class__.__name__} #{self.identifier}>"


class Boundary(ZBounded, ABC):
    """
    Can tell whether a Particle is inside of it or not.

    Making the total set of Boundary objects in an experiment closed is - for now - the job of the person implementing
    the experiment setup.
    """
    @abstractmethod
    def is_particle_inside(self, position: np.array, time: float) -> bool:
        """
        Tells whether the given position is inside this boundary at the given time.

        :param position: A spatial position to tell this for.
        :param time: The time to tell this for.
        :return: True if the particle is definitely inside this Boundary,
            False if it is not inside this Boundary (or if this is unknown).
        """
        pass


class Detector(ZBounded, ABC):
    """
    Can tell whether a Particle has hit it, and _if_ it has hit it, can also tell where this occurred.

    As opposed to real detectors like fluorescent screens, these Detectors must never have an effect on particles.
    This allows us to, for instance, make an image of a particle distribution at different places of the experiment,
    without terminating the path of any particle.

    Detector subclasses can - and should - freely choose how to internally do the calculation of the hit position.
    """
    def __init__(self, identifier: str):
        """
        The constructor for Detector.

        :param identifier: A unique identifier for this detector. Used to find which hit occurred on which detector.
        """
        if type(identifier) is not str:
            logging.warning(f'The passed detector identifier "{identifier}" is not a string, autoconverting to one...')
            identifier = str(identifier)
        self.identifier = identifier
        super().__init__()

    def try_to_detect(self, particle: Particle) -> bool:
        """
        Tries to detect a particle.

        If the particle is considered having hit the detector ("could be detected")
        based on its current state, then a DetectorHit is stored on the particle's `detector_hits` property
        and True is returned. Otherwise, this function has no effect and returns False.

        :param particle: A Particle instance.
        :return: True if the particle has reached this detector, False otherwise.
        """
        if self._has_particle_reached_detector(particle.identifier, particle.position):
            hit_position = self._hit_position(particle.phase_space_position)
            hit = ParticleDetectorHit(hit_position=hit_position, particle=particle)
            if self.identifier in particle.detector_hits:
                particle.detector_hits[self.identifier].append(hit)
            else:
                particle.detector_hits[self.identifier] = [hit]
            return True
        else:
            return False

    @abstractmethod
    def _has_particle_reached_detector(self, particle_identifier: int, position: np.array) -> bool:
        """Tells whether a Particle has reached this Detector. Must return True/False."""
        pass

    @abstractmethod
    def _hit_position(self, phase_space_position: np.array) -> Optional[np.array]:
        """
        Can return a hit position on this detector for a Particle, but might
        also return None (if the Particle isn't considered having reached this detector).

        :param phase_space_position: Position and velocity of the particle.
        :return: An (n,)-shaped numpy array describing the (spatial) hit position if there is a hit position
            to calculate, None otherwise. n is the number of spatial dimensions in the experiment.
        """
        pass


class Field(ZBounded, ABC):
    """
    A virtual acceleration field. Interacts with particles by exerting an acceleration on them, based on some
    collection of (local and current) properties of the particle and the field.
    """
    @abstractmethod
    def calculate_acceleration(self, particle: Particle, time: float) -> np.array:
        """
        Calculates an acceleration for one particle based on the particle's current properties and the current time.
        This acceleration will be integrated for in each time step and thus "applied" to the particle.

        :param particle: The Particle to calculate this Field's acceleration for.
        :param time: The time to calculate the acceleration for.
        :return: A (n,)-shaped numpy array describing the acceleration exerted on the particle. n is the number of
            spatial dimensions of the experiment.
        """
        pass


class Action(ABC):
    """
    An object to update a Particle's properties based on its current state. Can do arbitrary calculations to determine
    the values the properties should be set to.

    .. warning::
        If the properties that a .__call__() call changes include the .position attribute, the __call__ method
        MUST return True, and it SHOULD return False otherwise. If your action doesn't adhere to
        this, expect these updated positions to have no effect, since the integrator will overwrite them.

    Should be used for all calculations of additional quantities beyond calculating an acceleration, so,
    beyond what a Field returns. For example, TrackTrajectory is a simple Action that just takes the
    position and appends it to the `trajectory` property on the particle.

    The `__call__` method is guaranteed by `Experiment` to be called exactly once per time step (as long as the
    particle is still being simulated).

    If you need access to properties of fields, detectors, etc., override the __init__ method and store a reference
    to each relevant object at construction of the Action instance.

    .. note::
        Taking care to avoid data clashes (different actions writing on the same property on a particle,
        overwriting each others' results) is up to the user / implementer. This cannot be avoided in a general way. If
        you know that another Action instance is writing to a specific particle property, don't overwrite it unless you
        know the order of execution of the Actions and are doing this completely intentionally.
    """
    @abstractmethod
    def __call__(self, particle: Particle, time: float) -> bool:
        """
        Does some thing with a Particle instance. Must return True if the particle's ``phase_space_position`` is
        changed, or if its ``position`` or ``velocity`` is changed (since these two properties are directly derived
        from ``phase_space_position``). Note that a "change" entails overwriting the whole property as well as changing
        only parts of it.

        :param particle: The Particle instance.
        :param time: The current time.
        :return: True if the phase space position of a particle was changed in any way. False otherwise.
        """
        pass


class Device(ZBounded, ConfigSubscriber, ABC):
    """
    A combination of Fields, Actions and a Boundary. Used to model real-world devices in an experiment setup.
    """

    def __init__(self, fields: List[Field], boundary: Boundary, actions: Union[List[Action]] = None):
        """
        Constructor for Device.

        :param fields: The Fields that are present in this device.
        :param boundary: The Boundary this device has.
        """
        self._fields: List[Field] = []
        self._auto_overwrite_calculate_acceleration()
        self._boundary = boundary
        for field in fields:
            self.add_field(field)
        self._actions: List[Action] = actions if actions is not None else []

        self._number_of_dimensions = None
        super().__init__()
        GlobalConfig().subscribe(self, ConfigKey.NUMBER_OF_DIMENSIONS)

    def config_change(self, key: ConfigKey, value: Any):
        if key is ConfigKey.NUMBER_OF_DIMENSIONS:
            self._number_of_dimensions = value

    def calculate_acceleration(self, particle: Particle, time: float) -> np.array:
        """
        Calculates and returns the acceleration that this Device exerts on a Particle ``particle`` at time
        ``time``.

        :param particle: The Particle to calculate the acceleration for.
        :param time: The time to calculate the acceleration for.
        :return: The acceleration vector.
        """
        raise Exception("This function should never have been called, "
                        "implementation should have been automatically swapped out or overridden by subclass!")

    def is_particle_inside(self, particle_position: np.array, time: float) -> bool:
        """
        Returns whether the passed Particle is inside this Device or not.

        :param particle_position: The particle's spatial position.
        :param time: The current time.
        :return: True if the Particle inside this Device, False if it is not (or if this is unknown).
        """
        return self._boundary.is_particle_inside(particle_position, time)

    def run_actions(self, particle: Particle, time: float) -> bool:
        """
        Runs all Actions stored on this Device, *if* the passed Particle is inside this device.
        Returns True if any of them changed the particle position.
        :param particle: The Particle that all of this Device's Actions will be run for.
        :param time: The current time.
        :return: True if any of the Actions returned True, False otherwise.
        """
        position_changed = False
        if self.is_particle_inside(particle.position, time):
            for action in self.actions:
                position_changed |= action(particle, time)
        return position_changed

    @property
    def z_boundary(self) -> Tuple[float, float]:
        """The Z boundary (z_min, z_max) of this device."""
        return self._boundary.z_boundary

    @property
    def fields(self) -> List[Field]:
        """The list of fields (:class:`Field` instances) that this Device contains."""
        return self._fields

    def add_field(self, field: Field) -> None:
        """
        Adds a field (:class:`Field` instance) to this Device.

        :param field: The field to add.
        """
        self._fields.append(field)
        self._auto_overwrite_calculate_acceleration()

    @property
    def actions(self) -> List[Action]:
        """The list of actions (:class:`Action` instances) that this Device contains."""
        return self._actions

    def add_action(self, action: Action) -> None:
        """
        Adds an action (:class:`Action` instance) to this Device.

        :param action: The action to add.
        """
        self._actions.append(action)

    @property
    def boundary(self):
        """
        The boundary (:class:`Boundary` instance) that this device has.

        :return: The current boundary of this device.
        """
        return self._boundary

    @boundary.setter
    def boundary(self, boundary: Boundary):
        self._boundary = boundary

    def _auto_overwrite_calculate_acceleration(self):
        """
        If calculate_acceleration hasn't been overridden in the concrete class, this method overwrites the
        implementation of calculate_acceleration with an optimal implementation based on the number of Fields set on
        this Device.
        """
        if self.calculate_acceleration.__func__ in \
            [Device.calculate_acceleration, Device._calculate_acceleration_no_fields,
             Device._calculate_acceleration_one_field, Device._calculate_acceleration_multiple_fields,
             Device._calculate_acceleration_many_fields]:
            if len(self.fields) == 0:
                self.calculate_acceleration = self._calculate_acceleration_no_fields
            elif len(self._fields) == 1:
                self.calculate_acceleration = self._calculate_acceleration_one_field
            elif len(self._fields) <= 10:
                self.calculate_acceleration = self._calculate_acceleration_multiple_fields
            else:
                self.calculate_acceleration = self._calculate_acceleration_many_fields

    # noinspection PyUnusedLocal
    def _calculate_acceleration_no_fields(self, particle: Particle, time: float) -> np.array:
        return np.zeros(self._number_of_dimensions)

    def _calculate_acceleration_one_field(self, particle: Particle, time: float) -> np.array:
        return self._fields[0].calculate_acceleration(particle, time)

    def _calculate_acceleration_multiple_fields(self, particle: Particle, time: float) -> np.array:
        # For n <= 10 fields, this is roughly 10 times faster than using:
        # np.sum([f.calculate_acceleration(particle,time) for f in self.fields], axis=0)
        acceleration = np.zeros(self._number_of_dimensions)
        for field in self._fields:
            acceleration += field.calculate_acceleration(particle, time)
        return acceleration

    def _calculate_acceleration_many_fields(self, particle: Particle, time: float) -> np.array:
        # For n > 10 fields, this implementation should typically be faster than the one for _multiple_fields.
        return np.sum([f.calculate_acceleration(particle, time) for f in self._fields], axis=0)


class ResultStorage(ABC):
    """
    An object to store the results of an experiment in some fashion. MUST implement store_results, and MAY implement
    convenience methods to read from the storage again (e.g. a method to get all particle trajectories from a file).
    """
    @abstractmethod
    def store_results(self, particles: List[Particle]):
        """
        Stores the results of an experiment (which are always a list of modified Particle instances).

        :param particles: The list of particles, each in the state of after running a simulation.
        """
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    @abstractmethod
    def get_dimensions(self) -> int:
        """
        The number of spatial dimensions that the stored simulation had.
        """
        pass

    @abstractmethod
    def get_properties(self) -> Optional[Iterable]:
        """
        The collection of constant particle properties, i.e., the values of their
        :meth:`Particle.constant_properties`.

        :return: See above; None if these weren't stored.
        """
        pass

    @abstractmethod
    def get_tracked_initial(self) -> Optional[Iterable]:
        """
        The collection of initial particle states, i.e., the initial values of their
        :meth:`Particle.tracked_properties`.

        :return: See above; None if these weren't stored.
        """
        pass

    @abstractmethod
    def get_tracked_final(self) -> Optional[Iterable]:
        """
        The collection of final particle states, i.e., the final values of their
        :meth:`Particle.tracked_properties`.

        :return: See above; None if these weren't stored.
        """
        pass

    @abstractmethod
    def get_trajectories(self) -> Optional[Iterable]:
        """
        The collection of all particles' trajectories.

        :return: See above; None if these weren't stored.
        """
        pass

    @abstractmethod
    def get_detectors(self) -> Optional[Dict[str, Iterable]]:
        """
        The collection of all detectors.

        :return: The detectors, in a dictionary mapping their identifying string to their collection of detected hits.

        .. note::
          Detectors that did not detect at least one particle are not necessarily part of the output of this method,
          and that if there were no detectors or no detector detected at least one particle, None is returned.
        """
        pass


class Setup(ABC):
    """
    A base class for classes that define experiment setups. Should be considered static, i.e. its method are all
    "@staticmethod"s, and no instances should (need to be) created.
    """
    @staticmethod
    @abstractmethod
    def get_parser() -> SetupArgumentParser:
        """
        Returns a parser for the arguments relevant to this specific setup.
        This parser will receive only the arguments that the main program (bin/cminject) did not recognise.

        :return: An argparse.ArgumentParser (or subclass) instance that can parse all the args relevant to this setup.
        """
        pass

    @staticmethod
    def validate_args(args: argparse.Namespace):
        """
        Validates the arguments. Useful for validation that needs to check multiple args at once, and not just one
        specific arg. Overriding this method is optional and only needs to be done if such validation is desired.

        :param args: The argparse.Namespace object that the parser constructed by get_parser() returned after being
            given all arguments that the main program did not recognise.
        :raises: argparse.ArgumentError if validation failed
        """
        pass

    @staticmethod
    @abstractmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> 'Experiment':
        """
        Constructs an Experiment instance from the arguments to the main program and this specific setup.

        :param main_args: The argparse.Namespace parsed by the main program. Contains args that are common
            to all setups, like the number of particles and timespan.
        :param args: The argparse.Namespace parsed by this class's parse_args method.
        :return: An Experiment instance that's ready to be ran.
        """
        pass


class Source(ABC):
    """
    A source of particles, to generate an initial position (and property) distribution.

    One can implement sources following different distributions and generating different particle types by subclassing
    this class and implementing generate_particles however they want.

    .. warning::
      Any implementation of Source should only generate exactly one type of particle, i.e. its generate_particles
      method should never return a list of particles where the (most specific) type of any two particles is different.
    """
    @abstractmethod
    def generate_particles(self, start_time: float = 0.0) -> List[Particle]:
        """
        Generates a list of particles. How this is done is entirely up to the subclass, by (this is mandatory!)
        implementing this method in some way.

        :return: A list of Particle instances.
        """
        pass


class ParticleDetectorHit(object):
    """
    A small data container class, to store required information about a detection event of a particle by a detector.

    Contains :attr:`hit_state`, which is a full description of the particle's state, with the 'position' part
    overwritten by the position the detector detected the particle at.
    """
    def __init__(self, hit_position: np.array, particle: Particle):
        """
        The constructor for a ParticleDetectorHit.

        :param hit_position: An (n,)-shaped np.array representing the hit position of the particle on the detector.
            n is the number of spatial dimensions and must match that of the particle.
        :param particle: The particle which was detected.
        """
        dim_pos = hit_position.size
        dim_par = particle.number_of_dimensions
        if dim_pos != dim_par:
            raise ValueError(
                f"Hit position dimensionality ({dim_pos}) does not match the number of "
                f"dimensions the particle is simulated in ({dim_par})!"
            )

        hit_state = particle.as_array('all')
        hit_state['position'] = hit_position
        self._hit_state = hit_state

    @property
    def hit_state(self) -> np.array:
        """
        The state of the particle in the moment where the detector detected the particle.
        The 'position' part of this state refers to the (approximate) position of detection as returned by the detector.

        :return: The position of the hit.
        """
        return self._hit_state

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
