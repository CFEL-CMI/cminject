from typing import Dict, Any

import numpy as np
from cminject.experiment_refactored.definitions.base import PropertyUpdater, Particle
from cminject.experiment_refactored.definitions.particles import ThermallyConductiveSphericalParticle
from cminject.experiment_refactored.fields.fluid_flow_field import FluidFlowField


class TrajectoryPropertyUpdater(PropertyUpdater):
    def update(self, particle: Particle, time: float) -> None:
        particle.trajectory.append(
            np.concatenate([[time], particle.position])
        )


class DensityMapPropertyUpdater(PropertyUpdater):
    def __init__(self,
                 results: Dict[str, Any],
                 extent: np.array,
                 delta: float):
        bins = np.ceil(extent / delta).astype('int64')
        shape = bins[:, 1] - bins[:, 0]
        results['density_map'] = np.zeros(shape, dtype='int64')

        self.results = results
        self.extent = extent
        self.shape = shape
        self.delta = delta

        from scipy.interpolate import interp1d
        interpolators = []
        for i, axis in enumerate(extent):
            x = axis
            y = np.array([0, shape[i]-1])
            interp = interp1d(x, y, kind='linear')
            interpolators.append(interp)
        self.interpolators = interpolators

    def update(self, particle: Particle, time: float) -> None:
        bin_ = np.round([
            self.interpolators[i](particle.position[i])
            for i in range(len(particle.position))
        ]).astype(int)
        try:
            self.results['density_map'][bin_[0], bin_[1], bin_[2]] += 1
        except IndexError:
            pass


class ParticleTemperaturePropertyUpdater(PropertyUpdater):
    def __init__(self, field: FluidFlowField, dt: float):
        self.field: FluidFlowField = field
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