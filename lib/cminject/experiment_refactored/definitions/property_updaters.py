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


class MirrorSymmetryPropertyUpdater(PropertyUpdater):
    """
    A property updater that ensures mirror symmetry around an axis, specified by a dimension  and position relative
    to the origin in that dimension.

    Useful for e.g. radially symmetric 2D simulations around a central axis, where particles are not supposed
    to be simulated with negative r.

    It does this by checking the position relative to the specified axis and -- if the particle position is less than
    the axis position -- flipping position around the axis and flipping the dimension's velocity component.
    """
    def __init__(self, dimension: int = 0, position: float = 0.0):
        """
        The constructor for MirrorSymmetryPropertyUpdater.
        :param dimension: The index of the dimension. 0 by default, i.e. the first dimension.
        :param position: The position of the mirror axis relative to the origin. 0.0 by default, i.e. the default
        axis goes through the origin.
        """
        self.dimension = dimension
        self.position = position
        self.number_of_dimensions = 2

    def update(self, particle: Particle, time: float) -> None:
        if particle.position[self.dimension] < self.position:
            # Flip the position around the mirror axis
            particle.position[self.dimension] = 2 * self.position - particle.position[self.dimension]
            # Flip the velocity
            particle.position[self.dimension + self.number_of_dimensions] *= -1

    def set_number_of_dimensions(self, number_of_dimensions: int):
        self.number_of_dimensions = number_of_dimensions


class ParticleTemperaturePropertyUpdater(PropertyUpdater):
    """
    A property updater to calculate particle temperatures based on two models for particles within a
    StokesFluidFlowField.

    WARNING: This is untested and most likely incorrect at the moment. We'll need to revisit the code in the flow field
    that does the calculations if this code becomes necessary to use.
    """
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
