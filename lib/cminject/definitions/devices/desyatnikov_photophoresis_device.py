from typing import List

import numpy as np
from scipy.constants import Boltzmann

from cminject.definitions import Device, Field, PropertyUpdater
from cminject.definitions.boundaries import CuboidBoundary
from cminject.definitions.fields.laser_fields import DesyatnikovPhotophoreticLaserField
from cminject.definitions.fields.function_field import FunctionField
from cminject.definitions.particles import ThermallyConductiveSphericalParticle


# TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TODO this and the property updater should not live here! They should either not exist, or fall within
# TODO the family of drag force fields (just without grid interpolation). I chose to put them here for now, to allow
# TODO my combined simulations to finish faster without having to implement something new. This isn't good design,
# TODO it's only for me to be able to finish the thesis.     - Simon
# TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def get_uniform_drag_force(viscosity, velocity, temperature, m_gas, pressure):
    def uniform_drag_force(particle, time):
        knudsen = viscosity / (pressure * particle.radius) * \
                  np.sqrt(np.pi * Boltzmann * temperature / (2 * m_gas))
        s = 1 + knudsen * (1.231 + 0.4695 * np.exp(-1.1783 / knudsen))
        relative_velocity = velocity - particle.velocity

        force = 6 * np.pi * viscosity * particle.radius * relative_velocity * s
        return force / particle.mass
    return uniform_drag_force


class UniformBrownianMotionPropertyUpdater(PropertyUpdater):
    def __init__(self, viscosity, temperature, m_gas, pressure, dt):
        self.viscosity = viscosity
        self.temperature = temperature
        self.m_gas = m_gas
        self.pressure = pressure
        self.dt = dt

    def update(self, particle: ThermallyConductiveSphericalParticle, time) -> bool:
        if self.pressure >= 0.0:
            knudsen = self.viscosity / (self.pressure * particle.radius) * \
                      np.sqrt(np.pi * Boltzmann * self.temperature / (2 * self.m_gas))
            s = 1 + knudsen * (1.231 + 0.4695 * np.exp(-1.1783 / knudsen))

            s0 = 216 * self.viscosity * Boltzmann * self.temperature /\
                 (np.pi**2 * (2 * particle.radius)**5 * particle.rho**2 * s)
            a = np.random.normal(0.0, 1.0, 2) * np.sqrt(np.pi * s0 / self.dt)

            position = particle.spatial_position + (0.5 * a * self.dt**2)
            velocity = particle.velocity + (a * self.dt)
            particle.position = np.concatenate([position, velocity])

        return True

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != 2:
            raise ValueError("Only available in n=2 dimensions!")


class DesyatnikovVortexLaserDevice(Device):
    def __init__(self, r_boundary, z_boundary,
                 gas_temperature, gas_viscosity, gas_thermal_conductivity, gas_density, gas_pressure, gas_mass,
                 beam_power, beam_waist_radius, drag_field_velocity=None, z_position=0.0):

        if beam_power:
            pp_field = DesyatnikovPhotophoreticLaserField(
                gas_viscosity=gas_viscosity, gas_temperature=gas_temperature,
                gas_thermal_conductivity=gas_thermal_conductivity, gas_density=gas_density,
                beam_power=beam_power, beam_waist_radius=beam_waist_radius, z_position=z_position
            )
            fields: List[Field] = [pp_field]
        else:
            fields = []

        if drag_field_velocity:
            # Assume air at 100 mbar
            drag_field = FunctionField(get_uniform_drag_force(gas_viscosity, np.array([0, drag_field_velocity]),
                                                              gas_temperature, gas_mass, gas_pressure))
            fields.append(drag_field)

        boundary: CuboidBoundary = CuboidBoundary(intervals=[r_boundary, z_boundary])
        super().__init__(fields=fields, boundary=boundary)