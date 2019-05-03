from typing import List
from multiprocessing import Pool, Manager
from functools import partial

import numpy as np
from scipy.integrate import ode

from .base_classes import Particle, Source, Device, Detector


def calculate_v_and_a(time: float, position_and_velocity: np.array,
                      particle: Particle, devices: List[Device]):
    total_acceleration = np.zeros(3, dtype=float)
    for device in devices:
        if device.field.needs_full_particle_data:
            total_acceleration += device.field.calculate_acceleration(position_and_velocity, particle)
        else:
            total_acceleration += device.field.calculate_acceleration(position_and_velocity)
    return np.concatenate((position_and_velocity[3:], total_acceleration))


def simulate_particle(particle: Particle, devices: List[Device],
                      t_start: float, t_end: float, dt: float,
                      track_trajectories: bool):
    integral = ode(calculate_v_and_a)
    integral.set_integrator('lsoda')  # TODO maybe use other integrators?
    integral.set_initial_value(np.concatenate([particle.position, particle.velocity]), t_start)
    integral.set_f_params(particle, devices)

    while integral.successful() and integral.t < t_end:
        if True:  # TODO check if particle is not yet lost (outside the boundary of the whole experiment)
            particle.position = integral.y[:3]
            particle.velocity = integral.y[3:]
            integral.integrate(integral.t + dt)
        else:
            break

    return particle


class Experiment:
    def __init__(self, devices: List[Device], sources: List[Source], detectors: List[Detector],
                 dt: float = 1.0e-5, t_start=0.0, t_end=1.8, track_trajectories=False):
        # Store references
        self.devices = devices
        self.sources = sources
        self.detectors = detectors

        # Store numerical data
        self.dt = dt
        self.t_start = t_start  # TODO needed?
        self.t_end = t_end
        self.track_trajectories = track_trajectories

        # Initialise list of particles from all sources
        self.particles = []
        for source in self.sources:
            source_particles = source.generate_particles()
            self.particles += source_particles

    def run(self):
        print("Running experiment...")
        pool = Pool()
        try:
            parallel_simulate = partial(
                simulate_particle,
                devices=self.devices,
                t_start=self.t_start,
                t_end=self.t_end,
                dt=self.dt,
                track_trajectories=self.track_trajectories
            )  # TODO further params
            result = pool.imap(parallel_simulate, self.particles, chunksize=32)
        except Exception as e:
            raise e
        finally:
            pool.close()
            pool.join()

        return list(result)
