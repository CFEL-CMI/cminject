import numpy as np
from cminject.experiment_refactored.definitions.base import PropertyUpdater, Particle
from cminject.experiment_refactored.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.definitions.fields.stokes_fluid_flow_field import StokesFluidFlowField


class TrajectoryPropertyUpdater(PropertyUpdater):
    def update(self, particle: Particle, time: float) -> None:
        particle.trajectory.append(
            np.concatenate([[time], particle.position])
        )

    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass


class RadialSymmetryPropertyUpdater(PropertyUpdater):
    def __init__(self, dimension: int = 0, position: float = 0.0):
        self.dimension = dimension
        self.position = position
        self.number_of_dimensions = 2

    def update(self, particle: Particle, time: float) -> None:
        if particle.position[self.dimension] < 0.0:
            particle.position[self.dimension] = abs(particle.position[self.dimension])
            particle.position[self.dimension + self.number_of_dimensions] *= -1

    def set_number_of_dimensions(self, number_of_dimensions: int):
        self.number_of_dimensions = number_of_dimensions


class ParticleTemperaturePropertyUpdater(PropertyUpdater):
    def __init__(self, field: StokesFluidFlowField, dt: float):
        self.field: StokesFluidFlowField = field
        self.dt: float = dt

    def update(self, particle: ThermallyConductiveSphericalParticle, time: float) -> None:
        if self.field.is_particle_inside(particle):
            t = self.field.calc_particle_temperature(particle, time)
            n_c, t_c = self.field.calc_particle_collisions(particle, self.dt, time)

            particle.collision_temperature = t_c
            particle.temperature = t

            if t <= 77.0 and not np.isfinite(particle.time_to_liquid_n):
                particle.time_to_liquid_n = time
            if t_c <= 77.0 and not np.isfinite(particle.collision_time_to_liquid_n):
                particle.collision_time_to_liquid_n = time

    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass
