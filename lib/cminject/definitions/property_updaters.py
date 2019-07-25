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

from typing import Tuple

import numpy as np

from scipy.constants import pi, Avogadro, Boltzmann, R

from cminject.definitions.base import PropertyUpdater, Particle
from cminject.definitions.particles import ThermallyConductiveSphericalParticle, SphericalParticle
from cminject.definitions.fields.stokes_fluid_flow_field import StokesFluidFlowField
from cminject.definitions.fields.molecular_flow_drag_field import MolecularFlowDragField


class TrajectoryPropertyUpdater(PropertyUpdater):
    """
    A simple property updater to append to the trajectory of the particle after each time step.
    """
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


class BrownianMotionPropertyUpdater(PropertyUpdater):
    """
    Models brownian motion based on the paper:
    A. Li, G. Ahmadi, Dispersion and deposition of spherical particles from point sources in a turbulent channel flow,
    Aerosol Sci. Techn. 16 (24) (1992) 209–226. doi:10.1080/02786829208959550
    """
    def set_number_of_dimensions(self, number_of_dimensions: int):
        self.number_of_dimensions = number_of_dimensions  # TODO validate against field dimensions

    def __init__(self, field: StokesFluidFlowField, dt: float):
        self.field = field
        self.dt = dt
        self.number_of_dimensions = None

    def update(self, particle: SphericalParticle, time: float) -> None:
        # TODO what is this model ?
        _, _, pressure = self.field.interpolate(particle, time)
        if pressure >= 0.0:
            cunningham = self.field.calc_slip_correction(pressure, particle)
            s0 = 216 * self.field.dynamic_viscosity * Boltzmann * self.field.temperature /\
                 (pi**2 * (2 * particle.radius)**5 * particle.rho**2 * cunningham)
            a = np.random.normal(0.0, 1.0, self.number_of_dimensions) * np.sqrt(pi * s0 / self.dt)

            position = particle.spatial_position + (0.5 * a * self.dt**2)
            velocity = particle.velocity + (a * self.dt)
            particle.position = np.concatenate([position, velocity])


class BrownianMotionMolecularFlowPropertyUpdater(PropertyUpdater):
    """
    Models brownian motion based on the fluctuation-dissipation-theorem 
    and a numerical representation of the delta function
    """
    def set_number_of_dimensions(self, number_of_dimensions: int):
        self.number_of_dimensions = number_of_dimensions  # TODO validate against field dimensions

    def __init__(self, field: MolecularFlowDragField, dt: float):
        self.field = field
        self.dt = dt
        self.number_of_dimensions = None

    def update(self, particle: SphericalParticle, time: float) -> None:
        # TODO what is this model ?
        _, _, pressure = self.field.interpolate(particle, time)
        if pressure >= 0.0:
            h = self.field.m_gas / (2 * Boltzmann * self.field.temperature)
            s0 = (16/3 + 3*pi/5) * np.sqrt(pi/h) * pressure * self.field.m_gas * particle.radius**2
            if self.number_of_dimensions == 2:  # TODO this is just for radial symmetry
                a = np.random.normal(0.0, 1.0, 3) * np.sqrt(s0 / self.dt) / particle.mass
                particle.position[0] = np.sqrt((particle.position[0] + (0.5 * a[0] * self.dt**2))**2 +
                                               (0.5 * a[1] * self.dt**2)**2)
                particle.position[1] = particle.position[1] + (0.5 * a[2] * self.dt**2)
                particle.position[2] = np.sqrt((particle.position[2] + (a[0] * self.dt))**2 +
                                               (a[1] * self.dt)**2)
                particle.position[3] = particle.position[3] + (a[2] * self.dt)
            else:
                a = np.random.normal(0.0, 1.0, self.number_of_dimensions) * np.sqrt(s0 / self.dt) / particle.mass                
                position = particle.spatial_position + (0.5 * a * self.dt**2)
                velocity = particle.velocity + (a * self.dt)
                particle.position = np.concatenate([position, velocity])


class ParticleTemperaturePropertyUpdater(PropertyUpdater):
    """
    A property updater to calculate particle temperatures based on two models for particles
    within a StokesFluidFlowField.

    WARNING: This is untested and most likely incorrect at the moment. We'll need to revisit the code in the flow field
    that does the calculations if this code becomes necessary to use.
    """
    def __init__(self, field: StokesFluidFlowField, dt: float):
        self.field: StokesFluidFlowField = field
        self.dt: float = dt

    def update(self, particle: ThermallyConductiveSphericalParticle, time: float) -> None:
        if self.field.is_particle_inside(particle):
            t = self.calc_particle_temperature(particle, time)
            n_c, t_c = self.calc_particle_collisions(particle, self.dt, time)

            particle.collision_temperature = t_c
            particle.temperature = t

            if t <= 77.0 and not np.isfinite(particle.time_to_liquid_n):
                particle.time_to_liquid_n = time
            if t_c <= 77.0 and not np.isfinite(particle.collision_time_to_liquid_n):
                particle.collision_time_to_liquid_n = time

    def calc_thermal_conductivity(self, particle: Particle, time: float) -> float:
        """
        Calculates the thermal conductivity based on the ideal gas law
        http://hyperphysics.phy-astr.gsu.edu/hbase/thermo/thercond.html

        :param particle: A Particle
        :param time: The current time
        :return: The thermal conductivity
        """
        v: float = self.field.interpolate(particle, time)[1]
        mean_free_path = self.calc_mean_free_path()
        k = self.field.pressure * v * mean_free_path * self.field.molar_heat_capacity / \
            (3 * self.field.molar_mass * self.field.specific_gas_constant * self.field.temperature)
        return k

    def calc_prandtl(self, particle: Particle, time: float):
        """
        Calculates the Prandtl number.

        :param particle: A Particle
        :param time: The current time
        :return: The Prandtl number.
        """
        return self.field.specific_heat_capacity * self.field.dynamic_viscosity\
               / self.calc_thermal_conductivity(particle, time)

    def calc_reynolds(self, particle: Particle, d: float, time: float) -> float:
        """
        Calculates the Reynolds number.

        :param d: A characteristic length (in meters)
        :param particle: A Particle
        :param time: The current time.
        :return: The Reynolds number.
        """
        rho = self.field.pressure * self.field.specific_gas_constant * self.field.temperature
        v_abs = self.field.interpolate(particle, time)[1]
        return rho * v_abs * d / self.field.dynamic_viscosity

    def calc_nusselt(self, reynolds: float, prandtl: float, mu_s: float) -> float:
        """
        Calculates the Nusselt number based on the Reynolds and Prandtl number.

        :param mu_s: The dynamic viscosity at the surface temperature of nanoparticle
        :param reynolds: The Reynolds number.
        :param prandtl: The Prandtl number.
        :return: The Nusselt number
        """
        re, pr = reynolds, prandtl  # to keep the equation shorter
        nu = 2 + (0.4 * re**0.5 + 0.06 * re**0.67) * pr**0.4 * (self.field.dynamic_viscosity / mu_s)**0.25
        return nu

    def calc_mach(self, particle: Particle, time: float):
        """
        Calculates the Mach number of the fluid at a Particle's position

        :param particle: A Particle
        :param time: The current time
        :return: The Mach number
        """
        v: float = self.field.interpolate(particle, time)[1]
        return v / self.field.speed_of_sound

    def calc_convective_heat_transfer_coefficient(self, particle: Particle, nu: float, d: float, time: float):
        return nu / d * self.calc_thermal_conductivity(particle, time)

    def calc_particle_temperature(self, particle: ThermallyConductiveSphericalParticle, time: float):
        """
        Calculates a particle's "temperature" within this fluid at a given time.

        :param particle: A particle
        :param time: A time
        :return: The particle's temperature
        """
        d = 2 * particle.radius
        mu_s = 1.96e-5
        re = self.calc_reynolds(particle, d, time=time)
        pr = self.calc_prandtl(particle, time=time)
        m = self.calc_mach(particle, time=time)

        nu_0 = self.calc_nusselt(re, pr, mu_s)
        nu = nu_0 / (1 + 3.42 * (m * nu_0 / (re * pr)))

        h = self.calc_convective_heat_transfer_coefficient(particle, nu, d, time)
        c = particle.thermal_conductivity * particle.mass
        r = h * (4 * pi * particle.radius**2) / c  # TODO move surface area to particle class? method or field?
        return self.field.temperature + (particle.temperature - self.field.temperature) * np.exp(-r * time)

    def calc_particle_collisions(self,
                                 particle: ThermallyConductiveSphericalParticle,
                                 dt: float,
                                 time: float) -> Tuple[float, float]:
        """
        Estimates the number of collisions of a particle within a time frame based on Kinetic Gas Theory

        :param particle: A Particle
        :param dt: A time delta
        :param time: The current time
        :return: A tuple of (the number of collisions, the new temperature based on them)
        """
        v = np.sqrt(8 * self.field.temperature * Boltzmann / (self.m_gas * pi))
        v += self.field.interpolate(particle, time)[1]
        n = self.field.pressure * Avogadro / (R * self.field.temperature)
        c = 0.5 * n * v * dt * (4 * pi * particle.radius**2)
        k = (particle.mass + self.field.m_gas)**2 / (2 * particle.mass * self.field.m_gas)
        # TODO in the original code, dT was calculated but unused?
        collision_temperature = self.field.temperature + \
                                (particle.collision_temperature - self.field.temperature) * np.exp(-c / k)
        return c, collision_temperature

    def calc_mean_free_path(self):
        k = 1.38065e-23
        return self.field.dynamic_viscosity * np.sqrt(pi * k * self.field.temperature / (2 * self.field.m_gas))\
               / self.field.pressure

    def set_number_of_dimensions(self, number_of_dimensions: int):
        pass  # This class doesn't care about dimensionality at all


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""