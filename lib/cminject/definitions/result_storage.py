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
import os
from typing import List, Tuple

import h5py
import numpy as np

from .base import ResultStorage, Particle


class HDF5ResultStorage(ResultStorage):
    """
    Stores the results of an Experiment in an HDF5 file. The following keys will be available:

    - particles/{i}:

      - initial_position: The initial phase space position of the particle.
      - final_position: The final phase space position of the particle.
      - trajectory: The trajectory of the particle, which is an (n, d)-shaped array with an attached `description`
          string attribute containing names/descriptions for each column delimited by commas. n is the number of points
          along the trajectory, and d is the number of physical quantities (position, velocity, temperature, etc.)
          stored at each point.

    - detector_hits/{j}: All hits on the detector j. What's actually stored here is based on the array of the
        `full_properties` field of each `ParticleDetectorHit`, so depends on the implementation used. The default
        implementation for `ParticleDetectorHit` stores the concatenated array of [position, properties] of the particle
        at each hit. A `description` field should be attached similar to the one in `particles/i/trajectory`.

    The root node will also have a `dimensions` attribute attached, which is an integer denoting the number of spatial
    dimensions used in the simulation.

    .. note::
        For efficient and predictable storage, it's assumed that all detector hits on each detector return
        exactly the same shape of data, and the same property description. This should be the case in any sensible
        implementation of a detector anyways. If you do find you need to deviate from this, store an array of the
        maximum necessary length for each hit and fill out values you can't have for a specific hit with `np.nan`.
    """
    def __init__(self, filename: str, mode: str = 'r'):
        """
        Constructs a HDF5ResultStorage object.

        :param filename: The filename to store in/read from.
        :param mode: The mode to open the file with. 'r' by default to avoid overwriting data when using the convenience
            methods to read data from a written file, so if an instance should be used to actually store data,
            'w' or similar must be passed explicitly (refer to `h5py.File` docs for available modes).
        """
        self.filename = filename
        self.mode = mode
        # Verify that the output file doesn't already exist, let's avoid overwriting existing data
        if self.mode != 'r' and os.path.isfile(filename):
            raise ValueError("Output file already exists! Please delete or move the existing file.")

    def store_results(self, particles: List[Particle]) -> None:
        """
        Stores experiment results into an HDF5 file as described in the class docstring.

        :param particles: The list of particles as returned as the result from a finished Experiment run.
        :return: Nothing.
        """
        with h5py.File(self.filename, self.mode) as h5f:
            dimensions = len(particles[0].spatial_position)
            h5f.attrs['dimensions'] = dimensions
            detector_hits = {}

            for particle in particles:
                h5f[f'particles/{particle.identifier}/initial_position'] = np.array(particle.initial_position)
                h5f[f'particles/{particle.identifier}/final_position'] = np.array(particle.position)

                if particle.trajectory:
                    h5f[f'particles/{particle.identifier}/trajectory'] = np.array(particle.trajectory)
                    h5f[f'particles/{particle.identifier}/trajectory'].attrs['description'] = \
                        ['t'] + particle.position_description
                if particle.detector_hits:
                    for detector_id, hits in particle.detector_hits.items():
                        if detector_id not in detector_hits:
                            detector_hits[detector_id] = []
                        detector_hits[detector_id] += hits

            for detector_id, hits in detector_hits.items():
                if hits:
                    h5f[f'detector_hits/{detector_id}'] = np.array([hit.full_properties for hit in hits])
                    h5f[f'detector_hits/{detector_id}'].attrs['description'] = hits[0].full_properties_description

    def get_dimensions(self) -> int:
        """
        Gets the number of spatial dimensions stored on the HDF5 file. Useful for properly separating out e.g. positions
        and velocities from the arrays returned by get_trajectories and get_positions.

        :return: The number of spatial dimensions.
        """
        with h5py.File(self.filename, self.mode) as h5f:
            return int(h5f.attrs['dimensions'])

    def get_trajectories(self) -> List[np.array]:
        """
        Gets all stored trajectories as a list of np.arrays, each array being (d, n_i)-shaped, where d is the number
        of measured quantities of the stored results, and n_i is the number of data points along the trajectory
        of the i-th particle.

        :return: A list of np.arrays as described above.
        """
        trajectories = []

        with h5py.File(self.filename, self.mode) as h5f:
            for particle_id in h5f['particles']:
                if 'trajectory' in h5f[f'particles/{particle_id}'].keys():
                    trajectory = h5f[f'particles/{particle_id}/trajectory'][:].transpose()
                    trajectories.append(trajectory)

        return trajectories

    def get_initial_and_final_positions(self) -> Tuple[np.array, np.array]:
        """
        Gets all stored initial and final particle positions.

        :return: A 2-tuple containing:

            - One (d, m)-shaped np.array containing the initial positions
            - Another (d, m)-shaped np.array containing the final positions

            Here, d is the number of stored quantities for each particle, m is the number of stored particles.
        """
        with h5py.File(self.filename, self.mode) as h5f:
            particles = [h5f['particles'][k] for k in h5f['particles']]
            initial_positions = np.array([p['initial_position'] for p in particles]).transpose()
            final_positions = np.array([p['final_position'] for p in particles]).transpose()
            return initial_positions, final_positions

    def get_detectors(self) -> List[List[np.array]]:
        """
        Gets all stored detectors with their hits.

        :return: A list of np.arrays, each np.array corresponding to one detector and containing all hits as a
            (d, m_i)-shaped array.
        """
        detectors = []
        with h5py.File(self.filename, self.mode) as h5f:
            if 'detector_hits' in h5f:
                for detector_id in h5f['detector_hits']:
                    hits = h5f['detector_hits'][detector_id][:].transpose()
                    detectors.append(hits)
        return detectors


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End: