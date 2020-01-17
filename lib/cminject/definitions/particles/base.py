from abc import ABC, abstractmethod
from typing import Any, List, Dict

import numpy as np

from cminject.definitions.base import NDimensional


class Particle(NDimensional, ABC):
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
                 position: np.array, *args, **kwargs):
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
        self.detector_hits: Dict[int, List['ParticleDetectorHit']] = {}  # TODO naming is nonagnostic about detectors...
        self.mass: float = 0.0
        self.time_of_flight: float = start_time
        self.number_of_dimensions = None

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
        A list of strings matching the .properties attribute in length, describing each value in the array
        in some manner (most likely using standard physical abbreviations like rho, phi, T, k, ...).

        :return: See above.
        """
        pass

    @property
    def position_description(self) -> List[str]:
        """
        A list of strings matching the .position attribute in length, describing each value in the array
        in some manner (most likely using standard physical abbreviations like x, y, z, vx, vy, vz, ...).

        The default implementation returns:

        - ["x", "y", "z", "vx", "vy", "vz"] for 3D simulations
        - ["r", "z", "vr", "vz"] for 2D simulations
        - ["z", "vz"] for 1D simulations

        If your particle implementation deviates from this, override the property.

        :return: See above.
        """
        if self.number_of_dimensions == 3:
            return ['x', 'y', 'z', 'vx', 'vy', 'vz']
        elif self.number_of_dimensions == 2:
            return ['r', 'z', 'vr', 'vz']
        elif self.number_of_dimensions == 1:
            return ['z', 'vz']
        else:
            raise AttributeError(
                f"The default implementation of position_description is undefined for {self.number_of_dimensions} "
                f"dimensions!"
            )

    @property
    def reached_any_detector(self) -> bool:
        """
        Whether this particle has ever reached any detector.
        """
        return not not self.detector_hits  # `not not` to convert to boolean, to not leak data here

    @property
    def spatial_position(self) -> np.array:
        """
        The purely spatial position of the particle.
        """
        return self.position[:self.number_of_dimensions]

    @property
    def velocity(self) -> np.array:
        """
        The velocity of the particle.
        """
        return self.position[self.number_of_dimensions:]

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if self.position.size != 2 * number_of_dimensions:
            raise ValueError(
                f"Particle {self} is incompatible with {number_of_dimensions} dimensions: "
                f"Phase space position is 2*{self.position.size / 2} dimensional."
            )
        else:
            self.number_of_dimensions = number_of_dimensions

    def __str__(self):
        return f"<{self.__class__.__name__} #{self.identifier}>"