from typing import List

import numpy as np
from cminject.experiment_refactored.definitions.base import Particle


class SphericalParticle(Particle):
    """
    A simple spherical particle that has a radius and a density (rho).
    """

    def __init__(self, *args, radius: float, rho: float, **kwargs):
        super().__init__(*args, **kwargs)
        self.radius = radius
        self.rho = rho
        self.mass = self.rho * 4 / 3 * np.pi * (self.radius ** 3)
        self.number_of_dimensions = None

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions in [1, 2, 3]:
            super().set_number_of_dimensions(number_of_dimensions)
        else:
            raise ValueError("SphericalParticles can only be simulated in 1-, 2-, and 3D space.")

    @property
    def properties(self) -> np.array:
        return np.array([
            self.identifier,
            self.time_of_flight,
            self.radius,
            self.mass
        ])

    @property
    def properties_description(self) -> List[str]:
        return ['ID', 't', 'r', 'm']


class ThermallyConductiveSphericalParticle(SphericalParticle):
    """
    Lacking a better name (so far), this class extends on the SphericalParticle with a thermal conductivity
    that is required by e.g. the photophoretic force and for thermal calculations.
    """
    def __init__(self, *args,
                 thermal_conductivity: float = 6.3,
                 temperature: float = 298.0,
                 refraction_index: float = 1.0,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.thermal_conductivity = thermal_conductivity
        self.temperature = temperature
        self.collision_temperature = temperature
        self.refraction_index = refraction_index

        self.time_to_liquid_n = float('inf')
        self.collision_time_to_liquid_n = float('inf')

    @property
    def properties(self) -> np.array:
        basic_phase = super().properties
        return np.concatenate([
            basic_phase,
            [
                self.temperature,
                self.collision_temperature,
                self.time_to_liquid_n,
                self.collision_time_to_liquid_n,
            ]
        ])

    @property
    def properties_description(self) -> List[str]:
        return super().properties_description + [
            'T', 'T_c', 't_Ln', 't_Ln_c'
        ]
