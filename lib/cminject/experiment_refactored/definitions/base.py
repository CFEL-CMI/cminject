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

from abc import ABC, abstractmethod
from typing import Tuple, List, Any, Optional, Dict

import numpy as np

"""
An interval with no (with impossible) extent. Should be used to describe the boundary of an object
along an axis that has no extent along this axis.
"""
empty_interval = (float('inf'), float('-inf'))

"An interval with infinite extent."
infinite_interval = (float('-inf'), float('inf'))


class ParticleDetectorHit(object):
    """
    A small data container class, to store useful information about a hit of a particle on a detector.

    Contains :attr:`full_properties`, which is a full description of the particle's state starting with the hit position
    on the detector as a numpy array, and :attr:`full_properties_description`, which is a string list matching
    :attr:`full_properties` in size and describing each scalar stored in the array.
    """
    def __init__(self, hit_position: np.array, particle: 'Particle'):
        self.full_properties = np.concatenate([hit_position, particle.properties])
        self.full_properties_description = ['x_d', 'y_d', 'z_d'] + particle.properties_description

    @property
    def hit_position(self) -> np.array:
        """
        The position where the detector detected the particle hit.
        :return: The position of the hit.
        """
        return self.full_properties[:3]


class Particle(ABC):
    """
    Describes a particle whose trajectory we want to simulate.
    It is first and foremost a data container, and it and its subclasses should be written and used as such.

    It can be read by any part of the code, but should only be written to by instances of Experiment
    (or subclasses thereof).

    This class declares a few basic properties within __init__ that we expect every particle
    to have:
      - identifier (a unique id),
      - lost (a flag storing whether the particle is considered lost)
      - detector_hits (a dict mapping detector identifiers to hit lists)
      - position (the particle's position)
      - initial_position (the particle's initial position)
      - velocity (the particle's velocity)
      - mass (the particle's mass)
      - trajectory (a list describing points in the particle's path)
    """
    def __init__(self, identifier: Any, start_time: float,
                 position: np.array):
        """
        The constructor for Particle.
        :param identifier: The unique identifier the particle should have.
        :param start_time: The time the particle starts being simulated at.
        :param position: The phase space position (n-D position and velocity) the particle starts at.
               A (2*n,) dimensional vector, where n is the number of spatial dimensions.
        """
        self.identifier: int = identifier
        self.lost: bool = False
        self.position: np.array = np.copy(position)
        self.initial_position: np.array = np.copy(position)
        self.trajectory: List[np.array] = []
        self.detector_hits: Dict[int, List['ParticleDetectorHit']] = {}
        self.mass: float = self.calculate_mass()
        self.time_of_flight: float = start_time

    @abstractmethod
    def calculate_mass(self) -> float:
        """
        An abstract method that subclasses must implement to calculate the particle's mass.
        Expected to be called only once, at construction.
        :return: A floating-point value describing the mass in kg.
        """
        pass

    @property
    @abstractmethod
    def properties(self) -> np.array:
        """
        An abstract property that subclasses must implement to return the current properties of this particle:
        A (n,)-dimensional numpy array of all properties - beyond the usual phase space position (position+velocity) -
        that describe the particle's current state in a way that is useful to the problem domain.

        The result of this will be calculated and stored on each detector hit automatically,
        with the phase space position (particle.position) prepended to the front of the array.

        The size of the returned array MUST match the length of the list returned by `properties_description`.
        :return: A numpy array describing the particle's current phase.
        """
        pass

    @property
    @abstractmethod
    def properties_description(self) -> List[str]:
        """
        Returns a list of strings matching what .properties returns, describing each value in the array
        in some manner (most likely using standard physical abbreviations like rho, phi, T, k, ...).

        The size of the returned list must match the size of the numpy array returned by `properties`.
        :return: See above.
        """
        pass

    @property
    def reached_any_detector(self):
        return not not self.detector_hits  # `not not` to convert to boolean, to not leak data here

    def __str__(self):
        return f"<{self.__class__.__name__} #{self.identifier}>"


class Source(ABC):
    """
    A source of particles, to generate an initial phase space.

    One can implement sources following different distributions and generating different particle types by subclassing
    this class and implementing generate_particles however they want.
    """
    def __init__(self, position: np.array):
        """
        The constructor for Source.
        :param position: The x/y/z position where the source is placed.
        """
        self.position = position

    @abstractmethod
    def generate_particles(self, start_time: float = 0.0) -> List[Particle]:
        """
        Generates a list of particles. How this is done is entirely up to the subclass.
        The only important thing is adhering to the specification here and having this method always return
        a list of particles.
        :return: A list of Particle instances.
        """
        pass


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


class Detector(ZBounded, ABC):
    """
    Can tell whether a Particle has hit it, and _if_ it has hit it, can also tell where this occurred.

    As opposed to real detectors like fluorescent screens, these Detectors must never have an effect on particles.
    This allows us to, for instance, make an image of a particle distribution at different places of the experiment,
    without terminating the path of any particle.

    Detector subclasses can - and should - freely choose how to internally do the calculation of the hit position.
    """
    def __init__(self, identifier: int):
        """
        The constructor for Detector.
        :param identifier: A unique identifier for this detector. Used to find which hit occurred on which detector.
        """
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
        if self._has_particle_reached_detector(particle):
            hit_position = self._hit_position(particle)
            hit = ParticleDetectorHit(hit_position=hit_position, particle=particle)
            if self.identifier in particle.detector_hits:
                particle.detector_hits[self.identifier].append(hit)
            else:
                particle.detector_hits[self.identifier] = [hit]
            return True
        else:
            return False

    @abstractmethod
    def _has_particle_reached_detector(self, particle: Particle) -> bool:
        """Tells whether a Particle has reached this Detector. Must return True/False."""
        pass

    @abstractmethod
    def _hit_position(self, particle: Particle) -> Optional[np.array]:
        """
        Can return a hit position on this detector for a Particle, but might
        also return None (if the Particle isn't considered having reached this detector).
        :param particle: A Particle instance.
        :return: A (3,)-shaped numpy array describing the hit position if there is a hit position to calculate,
         None otherwise.
        """
        pass


class Boundary(ZBounded, ABC):
    """
    Can tell whether a Particle is inside of it or not.
    Making the total set of Boundary objects in an experiment closed is - for now - the job of the person implementing
    the experiment setup.
    """
    @abstractmethod
    def is_particle_inside(self, particle: Particle) -> bool:
        """
        Tells whether the passed Particle is inside of this Boundary or not.
        :param particle: The Particle to tell this for.
        :return: True if the particle is definitely inside this Boundary,
        False if it is not inside this Boundary or if this cannot be known.
        """
        pass


class Field(ZBounded, ABC):
    """
    A Field interacting with Particle objects.
    """
    @abstractmethod
    def calculate_acceleration(self,
                               particle: Particle,
                               time: float) -> np.array:
        """
        Calculates an acceleration for one particle based on the particle's current properties and the current time.
        This acceleration will be integrated for in each time step and thus "applied" to the particle.
        :param particle: The Particle to calculate this Field's acceleration for.
        :param time: The time to calculate the acceleration for.
        :return: A (3,)-shaped numpy array describing the x/y/z acceleration exerted on the particle.
        """
        pass


class Device(ZBounded, ABC):
    """
    A combination of one Field instance and one Boundary instance.
    Used to model real-world devices in an experiment setup.

    As the correspondence to Field and Boundary is one-to-one, complex kinds of Fields and Boundaries should be
    implemented internally in a subclass of these classes.
    """
    def __init__(self, field: Field, boundary: Boundary):
        """
        Constructor for Device.
        :param field: The Field this device has.
        :param boundary: The Boundary this device has.
        """
        self.field = field
        self.boundary = boundary
        super().__init__()

    def _calculate_z_boundary(self) -> Tuple[float, float]:
        """
        Returns the Z boundary of this device. Defaults to returning a Z boundary encompassing both
        the device's Z boundary and the field's Z boundary, but should be overridden if a different
        Z boundary is required.
        :return: A (z_min, z_max) tuple as defined in ZBoundedMixin.
        """
        field_boundary = self.field.z_boundary
        boundary = self.boundary.z_boundary
        return min(field_boundary[0], boundary[0]), max(field_boundary[1], boundary[1])

    def is_particle_inside(self, particle: Particle) -> bool:
        """
        Returns whether the passed Particle is inside this Boundary or not.
        Defaults to returning what the Boundary's method with the same name returns, but should be overridden
        if a more complex decision is required.
        :param particle: A Particle instance.
        :return: True if the Particle inside this device, False if it is not or if this is unknown.
        """
        return self.boundary.is_particle_inside(particle)


class PropertyUpdater(ABC):
    """
    An object to update a Particle's properties based on its current state. Can do arbitrary calculations to determine
    the values the properties should be set to.

    Should be used for all calculations of additional quantities beyond calculating an acceleration, so,
    beyond what a Field returns. For example, TrajectoryPropertyUpdater is a simple PropertyUpdater that just takes the
    position and appends it to the `trajectory` property on the particle.

    The :method:`update` method is guaranteed by `Experiment` to be called exactly once per time step (as long as the
    particle is still being simulated).

    If you need access to properties of fields, detectors, etc., override the __init__ method and store a reference
    to each relevant object at construction of the PropertyUpdater instance.

    NOTE:
        Taking care to avoid data clashes (different property updaters writing on the same property on a particle,
        overwriting each others' results) is up to the user / implementer. This cannot be avoided in a general way.
        If you know that another PropertyUpdater instance is writing to a specific particle property, don't write to it
        too, unless you know the order of execution of the PropertyUpdaters and are doing this completely intentionally.
    """
    @abstractmethod
    def update(self, particle: Particle, time: float) -> None:
        """
        Updates a property of some Particle instance in some way.
        :param particle: The Particle instance.
        :param time: The current time.
        :return: Nothing.
        """
        pass


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""