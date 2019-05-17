from typing import Tuple

import numpy as np

from cminject.experiment_refactored.base_classes import Field, Particle, ZBoundedMixin
from cminject.experiment_refactored.basic import infinite_interval, SphericalParticle
from scipy.interpolate import RegularGridInterpolator


class FluidFlowField(Field):
    def __init__(self, filename: str,
                 density: float, dynamic_viscosity: float,
                 temperature: float = 4.0,  # , pressure: float = 100.0, thermal_creep: float = 1.0,
                 # inflow_speed: float = 20.0, outflow_pressure: float = 0.0, k_n: float = 912.0,
                 # kinetic_d: float = 260.0e-12, conv: float = 0.00001,
                 # molar_mass: float = 0.004002602, molar_heat_capacity: float = 12.5,
                 # specific_gas_constant: float = 2077.0, specific_heat_capacity: float = 3116.0,
                 # m_gas: float = 6.6e-27, m_gas_mol: float = 0.004002602,
                 # speed_of_sound: float = 117.7,
                 scale_slip: float = 1.0):
        self.density = density
        self.dynamic_viscosity = dynamic_viscosity
        self.temperature = temperature
        self.scale_slip = scale_slip

        self.min_z, self.max_z = infinite_interval
        self.f_drag = None

        self.outside_particles = set()

        self.read_from_file(filename)

    def get_z_boundary(self) -> Tuple[float, float]:
        return self.min_z - 0.001, self.max_z + 0.001

    def is_particle_inside(self, particle: SphericalParticle) -> bool:
        return particle.identifier not in self.outside_particles

    def calculate_drag_force(self, particle: SphericalParticle) -> np.array:
        """
        This function calculates the drag force using Stokes' law for spherical particles in continuum
        """
        try:
            f_drag_result = self.f_drag((particle.position[0], particle.position[1], particle.position[2]))
            pressure = f_drag_result[3]
            if pressure <= 0:
                self.outside_particles.add(particle.identifier)
                return np.zeros(3)

            relative_velocity = f_drag_result[:3] - particle.velocity
            force_vector = 6 * np.pi * self.dynamic_viscosity * particle.radius * relative_velocity
            return force_vector / self.calculate_slip_correction(pressure=pressure, particle=particle)
        except ValueError:  # Can't calculate drag force, particle is outside the field boundaries
            self.outside_particles.add(particle.identifier)
            return np.zeros(3)

    def calculate_slip_correction(self, pressure: float, particle: SphericalParticle) -> float:
        """
        Calculates the slip correction factor with temperature corrections.
        The Sutherland constant for helium is 79.4 at reference temperature of 273.0 K.
        I took a reference pressure of 1 Pascal, in which the mean free path of helium is 0.01754.

        see J. Aerosol Sci. 1976 Vol. 7. pp 381-387 by Klaus Willeke
        """
        c_sutherland = 79.4
        ref_temperature = 273.0  # K

        factor = self.temperature / ref_temperature \
            * (1 + c_sutherland / ref_temperature) \
            / (1 + c_sutherland / self.temperature)
        knudsen = 0.01754 * 1 / (pressure * particle.radius) * factor

        s = 1 + knudsen * (1.246 + (0.42 * np.exp(-0.87 / knudsen)))
        return s * self.scale_slip

    def calculate_acceleration(self, particle: SphericalParticle, time: float) -> np.array:
        return self.calculate_drag_force(particle) / particle.mass

    @staticmethod
    def read_text(filename: str) -> Tuple[np.array, np.array, np.array, np.array, np.array, np.array, np.array]:
        print(f"Reading in {filename}...")
        f = open(filename)
        x, y, z = [], [], []
        vx, vy, vz = {}, {}, {}
        p = {}

        for i in f:
            token = i.split()
            if ('%' not in token) and ('NaN' not in token):
                x.append(float(token[0]))
                y.append(float(token[1]))
                z.append(float(token[2]))
                vx[(float(token[0]), float(token[1]), float(token[2]))] = float(token[3])
                vy[(float(token[0]), float(token[1]), float(token[2]))] = float(token[4])
                vz[(float(token[0]), float(token[1]), float(token[2]))] = float(token[5])
                p[(float(token[0]), float(token[1]), float(token[2]))] = float(token[6])

        x, y, z = [sorted(set(w)) for w in [x, y, z]]
        n_x, n_y, n_z = [len(w) for w in [x, y, z]]
        vx_out, vy_out, vz_out, p_out = [np.zeros((n_x, n_y, n_z)) for i in range(4)]

        for i in range(n_z):
            for j in range(n_y):
                for k in range(n_x):
                    try:
                        vx_out[k, j, i] = vx[(x[k], y[j], z[i])]
                        vy_out[k, j, i] = vy[(x[k], y[j], z[i])]
                        vz_out[k, j, i] = vz[(x[k], y[j], z[i])]
                        p_out[k, j, i] = p[(x[k], y[j], z[i])]
                    except KeyError:
                        pass

        print(f"Successfully read {filename}.")
        return x, y, z, vx_out, vy_out, vz_out, p_out

    def read_from_file(self, filename):
        x, y, z, v_x, v_y, v_z, p = self.read_text(filename)

        self.min_z = np.min(z)
        self.max_z = np.max(z)

        data_grid = np.zeros(v_x.shape + (4,))
        data_grid[:, :, :, 0] = v_x
        data_grid[:, :, :, 1] = v_y
        data_grid[:, :, :, 2] = v_z
        data_grid[:, :, :, 3] = p

        self.f_drag = RegularGridInterpolator((x, y, z), data_grid)
