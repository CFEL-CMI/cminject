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


class DetectorHit(object):
    def __init__(self, detector_identifier: int, hit_position: np.array, particle: "Particle"):
        self.particle = particle
        self.particle_phase = np.concatenate([particle.get_phase(), hit_position])
        self.particle_phase_description = particle.get_phase_description() + ['x_d', 'y_d', 'z_d']

        self.detector_identifier = detector_identifier  # TODO store ref to detector instead?
        self.hit_position = hit_position


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
        self.detector_hits: Dict[int, List[DetectorHit]] = {}
        self.position: np.array = np.copy(position)
        self.initial_position: np.array = np.copy(position)
        self.velocity: np.array = np.copy(velocity)
        self.initial_velocity: np.array = np.copy(velocity)
        self.trajectory: List[np.array] = []
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

    @abstractmethod
    def get_phase(self) -> np.array:
        """
        An abstract method that subclasses must implement to return the current phase of this particle, i.e. a
        numpy array of particle properties that describe the particle's current state in a way that is useful
        to the problem domain.

        This phase will be calculated and stored on each detector hit automatically (unless store_hit is overridden).

        The size of the returned array must match the length of the list returned by `get_phase_description`.
        :return: A numpy array describing the particle's current phase.
        """
        pass

    @abstractmethod
    def get_phase_description(self) -> List[str]:
        """
        Returns a list of strings matching what get_phase() returns, describing each value in the phase array
        in some manner (most likely using standard one-letter abbreviations like x,y,z,T,...).

        The size of the returned list must match the size of the array returned by `get_phase`.
        :return: See above.
        """
        pass

    def store_hit(self, detector_identifier: int, position: np.array) -> None:
        """
        Stores a hit on a detector within detector_hits.
        :param detector_identifier: The identifier of the detector. Should be unique, otherwise there'll be data clashes
        :param position: The position where the particle hit the detector.
        The particle is passed the position instead of calculating it because it's up to the detector to decide
        when and where a particle has hit it:
        Some detectors might do interpolation of the last position, and they might do it in different ways.
        :return: None
        """
        hit = DetectorHit(
            particle=self,
            detector_identifier=detector_identifier,
            hit_position=position
        )
        if detector_identifier in self.detector_hits:
            self.detector_hits[detector_identifier].append(hit)
        else:
            self.detector_hits[detector_identifier] = [hit]

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
    def __init__(self):
        self.z_boundary = self._calculate_z_boundary()

    @abstractmethod
    def _calculate_z_boundary(self) -> Tuple[float, float]:
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
        super().__init__()

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


class Calculator(ABC):
    """
    An object to do and store arbitrary calculations based on a Particle's properties and the current time.

    Should be used for all calculations of additional quantities beyond calculating an acceleration, so,
    beyond what a Field returns. For example, TrajectoryCalculator is a simple calculator that just takes the position
    and appends it to the `trajectory` property on the particle.

    The `calculate` method is guaranteed by `Experiment` to be called exactly once per time step, as long as the
    particle is still being simulated.

    Taking care to avoid data clashes (different calculators writing on the same property on a particle, overwriting
    each others' results) is up to the user / implementer.

    If you need access to properties of fields, detectors, etc., override the __init__ method and store a reference
    to each relevant object at construction of the Calculator instance.
    """
    @abstractmethod
    def calculate(self, particle: Particle, time: float) -> None:
        """
        Does a calculation. Should store its result on the particle.
        :param particle: The Particle in its current state.
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