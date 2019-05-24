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
from typing import Tuple, List, Any, Optional, Dict, Set

import numpy as np


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
                 position: np.array, velocity: np.array):
        """
        The constructor for Particle.
        :param identifier: The unique identifier the particle should have.
        :param start_time: The time the particle starts at.
        :param position: The position the particle starts at.
        :param velocity: The velocity the particle starts with.
        """
        self.identifier: int = identifier
        self.lost: bool = False
        self.detector_hits: Dict[int, List[np.array]] = {}
        self.position: np.array = np.copy(position)
        self.initial_position: np.array = np.copy(position)
        self.velocity: np.array = np.copy(velocity)
        self.initial_velocity: np.array = np.copy(velocity)
        self.trajectory: List[np.array] = []
        self.mass: float = self.calculate_mass()
        self.time: float = start_time

    @abstractmethod
    def calculate_mass(self) -> float:
        """
        An abstract method that subclasses must implement to calculate the particle's mass.
        Expected to be called only once, at construction.
        :return: A floating-point value describing the mass in kg.
        """
        pass

    def store_hit(self, identifier: int, position: np.array) -> None:
        """
        Stores a hit on a detector within detector_hits.
        :param identifier: The identifier of the detector. Should be unique, otherwise there'll be data clashes.
        :param position: The position where the particle hit the detector.
        The particle is passed the position instead of calculating it because it's up to the detector to decide
        when and where a particle has hit it:
        Some detectors might do interpolation of the last position, and they might do it in different ways.
        :return: Nothing.
        """
        if identifier in self.detector_hits:
            self.detector_hits[identifier].append(position)
        else:
            self.detector_hits[identifier] = [position]

    @property
    def reached(self):
        return self.detector_hits

    def __str__(self):
        return f"<{self.__class__.__name__} #{self.identifier}>"


class Source(ABC):
    """
    A source of particles. One can implement sources following different distributions
    by subclassing this class and implementing generate_particles however they want.
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


class ZBoundedMixin(ABC):
    """
    A mixin for objects in an experiment setup that are bounded in the Z direction.
    Classes deriving from this mixin will have to implement the get_z_boundary method.

    Objects implementing this mixin should be used to determine the minimum/maximum extent of the entire experiment
    setup in the Z direction, and then allows an Experiment to make a fast preliminary calculation about whether it
    can consider a particle lost based on this.
    """
    @abstractmethod
    def get_z_boundary(self) -> Tuple[float, float]:
        """
        Returns the Z boundary of this Z-bounded object.
        :return: A tuple of floats, the first entry being z_min, the second being z_max.
        """
        pass


class Detector(ZBoundedMixin, ABC):
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

    @abstractmethod
    def has_reached_detector(self, particle: Particle) -> bool:
        """
        Tells whether a Particle has reached this Detector. Must return True/False.
        :param particle: A Particle instance.
        :return: True if the particle has reached this detector, False otherwise.
        """
        pass

    @abstractmethod
    def get_hit_position(self, particle: Particle) -> Optional[np.array]:
        """
        Can return a hit position on this detector for a Particle, but might
        also return None (if the Particle isn't considered having reached this detector).
        :param particle: A Particle instance.
        :return: A (3,)-shaped numpy array describing the hit position if there is a hit position to calculate,
         None otherwise.
        """
        pass


class Boundary(ZBoundedMixin, ABC):
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


class Field(ZBoundedMixin, ABC):
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


class Device(ZBoundedMixin, ABC):
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

    def get_z_boundary(self) -> Tuple[float, float]:
        """
        Returns the Z boundary of this device. Defaults to returning a Z boundary encompassing both
        the device's Z boundary and the field's Z boundary, but should be overridden if a different
        Z boundary is required.
        :return: A (z_min, z_max) tuple as defined in ZBoundedMixin.
        """
        field_boundary = self.field.get_z_boundary()
        boundary = self.boundary.get_z_boundary()
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


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""