import numpy as np
from scipy.constants import Boltzmann

from cminject.definitions.fields.fluid_flow import MolecularFlowDragForceField
from .base import PropertyUpdater


class ParticleTemperaturePropertyUpdater(PropertyUpdater):
    def __init__(self, field: MolecularFlowDragForceField, dt: float):
        self.field = field
        self.dt = dt
        self.number_of_dimensions = None

    def set_number_of_dimensions(self, number_of_dimensions: int):
        self.number_of_dimensions = number_of_dimensions

    def update(self, particle: ThermallyConductiveSphericalParticle, time: float) -> bool:
        pressure = self.field.interpolate(particle.spatial_position)[self.number_of_dimensions]
        h = self.field.m_gas / (2 * Boltzmann * self.field.temperature)
        h_ = self.field.m_gas / (2 * Boltzmann * particle.temperature)
        deltaE = 4 * pressure * np.sqrt(np.pi) * particle.radius**2 * (np.sqrt(h)/h_ - 1/np.sqrt(h))
        deltaT = deltaE / (particle.specific_heat * particle.mass)

        particle.temperature += deltaT * self.dt
        return False