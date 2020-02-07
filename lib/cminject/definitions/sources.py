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

from typing import List, Type, Dict, Union

import numpy as np
from cminject.definitions.base import Source
from cminject.definitions.particles import SphericalParticle
from cmiinject.definitions.particles import Molecule
import h5py as hp

Distribution = Union[Dict, float]  # An appropriate type for variable distribution descriptions
# I think The class VariableDistributionSource should be edited to accept other kind of particles.

#ParticleType = Union[SphericalParticle, Molecule]

class VariableDistributionSource(Source):
    """
    A Source for particles that allows variable distributions for each dimension of the initial phase space.
    A "Distribution" is a dictionary describing a specific kind of distribution, with a 'kind' key that must match
    one of the following distributions, and further keys describing parameters for the given distribution:

    - gaussian: A gaussian distribution with keys 'mu' and 'sigma'. (like np.random.normal)
    - linear: A deterministically linear distribution with keys 'min' and 'max' (like np.linspace)
    - uniform: A random uniform distribution with keys 'min' and 'max' (like np.random.uniform)
    - radial_gaussian: A 1D-projected 2D gaussian distribution. Useful for a radial dimension.
    - radial_linear: An approximation of a 1D-projected 2D deterministically linear distribution.
        Useful for a radial dimension.
    - radial_uniform: A 1D-projected 2D random uniform distribution. Useful for a radial dimension.
    """
    def __init__(self,
                 number_of_particles: int,
                 position: List[Distribution],
                 velocity: List[Distribution],
                 radius: Distribution,
                 rho: float,
                 seed=None,
                 randomly_rotate_around_z: bool = False,
                 subclass: Type[SphericalParticle] = SphericalParticle,
                 **subclass_kwargs):
        self.number_of_particles = number_of_particles
        self.position = position
        self.velocity = velocity
        self.radius = radius
        self.rho = rho
        self.seed = seed
        self.randomly_rotate_around_z = randomly_rotate_around_z
        self.subclass = subclass
        self.subclass_kwargs = subclass_kwargs
        self.number_of_dimensions = len(self.position)

        # Check that position and velocity descriptions are equal in length (dimensionality)
        if len(self.position) != len(self.velocity):
            raise ValueError(
                f"Position and velocity descriptions must match in size! "
                f"Given position description was {len(self.position)}-dimensional, "
                f"velocity description was {len(self.velocity)}-dimensional."
            )

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"Incompatible number of dimensions: {number_of_dimensions}, "
                f"the position description passed at construction was {self.number_of_dimensions}-dimensional."
            )

    def _generate(self, dist: Distribution):
        if type(dist) is float:
            return np.repeat(dist, self.number_of_particles)

        if dist['kind'] == 'gaussian':
            mu, sigma = dist['mu'], dist['sigma']
            return np.random.normal(mu, sigma, self.number_of_particles)
        elif dist['kind'] == 'linear':
            l, r = dist['min'], dist['max']
            return np.linspace(l, r, self.number_of_particles)
        elif dist['kind'] == 'radial_gaussian':
            mu, sigma = dist['mu'], dist['sigma']
            random_x = np.random.normal(mu, sigma, self.number_of_particles)
            random_y = np.random.normal(0, sigma, self.number_of_particles)
            return np.sqrt(random_x**2 + random_y**2)
        elif dist['kind'] == 'radial_linear':
            l, r = dist['min'], dist['max']
            lin = np.linspace(l, r, self.number_of_particles)
            h = 1/np.sqrt(np.abs(lin[np.where(lin != 0)]))
            c = np.cumsum(h)
            return l + (c*(r-l) / np.sum(h))
        elif dist['kind'] == 'uniform':
            l, r = dist['min'], dist['max']
            return np.random.uniform(l, r, self.number_of_particles)
        elif dist['kind'] == 'radial_uniform':
            l, r = dist['min'], dist['max']
            random_x = np.random.uniform(l, r, self.number_of_particles)
            random_y = np.random.uniform(l, r, self.number_of_particles)
            return np.sqrt(random_x**2 + random_y**2)
        else:
            raise ValueError(f"Unknown or unspecified kind in distribution specification: {dist}")

    @staticmethod
    def _rotate_around_z(positions):
        size = positions.shape[0]
        rotated = (positions[:, 0] + positions[:, 1] * 1j) * np.exp(np.random.random(size) * 2 * np.pi * 1j)
        return np.array([rotated.real, rotated.imag, positions[:, 2]]).transpose()

    def generate_particles(self, start_time: float = 0.0):
        position = np.array([self._generate(pdist) for pdist in self.position]).transpose()
        if self.randomly_rotate_around_z:
            position = self._rotate_around_z(position)

        velocity = np.array([self._generate(vdist) for vdist in self.velocity]).transpose()
        r = self._generate(self.radius)

        particles = []
        for i in range(self.number_of_particles):
            inst = self.subclass(
                identifier=i,
                position=np.concatenate([position[i], velocity[i]]),
                start_time=start_time,
                radius=r[i],
                rho=self.rho,
                **self.subclass_kwargs
            )
            particles.append(inst)
        return particles

class MolDistributionSource(VariableDistributionSource):
    """
    Similar properties to the parent class but also includes Boltzmann distribution.
    """
    def __init__(self,  number_of_particles: int,
                position: List[Distribution],
                velocity: List[Distribution],
                radius: Distribution = 0,
                rho: float = 0,
                seed=None,
                randomly_rotate_around_z: bool = False,
                subclass: subclass: Type[Molecule] = Molecule,
                temperature: float,
                jmax : int,
                mass: float,
                energy_filename: str,
                **subclass_kwargs)):

        super().__init__(number_of_particles, position, velocity, radius, rho, seed, randomly_rotate_around_z, subclass, **subclass_kwargs)
        self.energy_filename = energy_filename
        self.temperature = temperature
        self.mass = mass
        self.jmax = jmax

    def BoltzmannAssigner(self, paths, ZeroFieldEnergies):
        """
        I still need to consider:
        1- The case in which the zero energy is 0
        2- nuclear statistics
        3- Add a maximum zero-field energy to consider
        4- Add tolerance
        Note: I should check later in what units are energies stored.
        """
        boltz = 1.38065053e-23
        Populations = []
        sum = 0
        for path, energy in zip(paths, ZeroFieldEnergies):
            pop = math.exp(-energy/(self.temperature*boltz))
            sum = sum + pop
            Populations.extend([pop])
        probabs = np.array(Populations)*(1/sum)
        ParticlesAssigned = np.ceil(probabs*self.number_of_particles)

        return ParticlesAssigned

    def PathFinder(self):
        energy_filename = self.energy_filename
        paths = []
        ZeroFieldEnergies = []
        qstates = [] # to store the quantum state of each file we want
        def pather(name, node):
            """
            The script finds the paths to all energies corresponding to j<=jmax. This function works together with hp.visititems()
            It also stores the quantum number of each state.
            Input:
            name: name of all paths. Usually outputed by hp.visititems
            Output:
            paths: a list of strings with all paths to the desired energies
            """
            jmax = self.jmax
            if 'energ' in name:
                slash_indices = []
                dash_indices = []
                for i, c in enumerate(name):
                    if c == '/':
                        slash_indices.append(i)
                    elif c == '_':
                        dash_indices.append(i)

                J = int(name[1:indices[0]])

                if J<=jmax:
                    paths.append(energy_filename)
                    Ka = int(name[dash_indices[1]+1, slash_indices[1]])
                    Kc = int(name[dash_indices[2]+1, slash_indices[2]])
                    M = int(name[dash_indices[3]+1, slash_indices[3]])
                    Isomer = int(name[dash_indices[4]+1, slash_indices[4]])
                    d = [J, Ka, Kc, M, Isomer]
                    qstates.append(d)
        with hp.File(energy_filename, 'r') as gl:
            gl.visititems(pather)
            # storing the free-field energy of each state.
            for path in paths:
                # energies in cmi_stark are not stored as numpy arrays. Should edit this later.
                energy = gl.get(path)
                energy = np.asarray(energy)
                energy = energy[0] # cuz it's stored as VLarray
                energy_0 = energy[0] # We're assuming the first point corresponds to a zero field.
                ZeroFieldEnergies.append(energy_0)
        return paths, ZeroFieldEnergies, np.array(qstates)

    def generate_particles(self, start_time: float = 0.0):
        position = np.array([self._generate(pdist) for pdist in self.position]).transpose()
        if self.randomly_rotate_around_z:
            position = self._rotate_around_z(position)
        velocity = np.array([self._generate(vdist) for vdist in self.velocity]).transpose()
        particles = []
        # Get the paths, free-field energies and quantum states of quantum states that have j < jmax
        paths, ZeroFieldEnergies, qstates = self.PathFinder()
        # Assign weights to each quantum state
        ParticlesAssigned = self.BoltzmannAssigner(paths, ZeroFieldEnergies)
        m_qstates = np.repeat(qstates, ParticlesAssigned, axis = 0)
        for i, q_state in zip(range(self.number_of_particles), m_qstates):
            # I need to construct the q_n several times depending on ParticlesAssigned
            d = {'J': q_state[0], 'Ka': q_state[1], 'Kc': q_state[2], 'M': q_state[3], 'Isomer':q_state[4]}

            inst = self.subclass(
                identifier=i,
                position=np.concatenate([position[i], velocity[i]]),
                start_time=start_time,
                mass=self.mass,
                q_n=d,
                e_filename = energy_filename
                **self.subclass_kwargs
            )
            particles.append(inst)
        return particles

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
