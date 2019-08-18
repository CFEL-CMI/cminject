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
from typing import List, Union, Tuple

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

    def get_trajectories(self, return_velocities: bool = False)\
            -> Union[List[np.array], Tuple[List[np.array], List[np.array]]]:
        """
        Returns all trajectories, which will be either one list of np.arrays (if `return_velocities` is False)
        containing each trajectory as an np.array of recorded positions, or a tuple of two such lists, the first
        containing each trajectory's positions and the second containing each trajectory's velocities (if
        `return_velocities` is True). If two lists are returned, they will match each other in length.

        :param return_velocities: Whether to also return the velocities in a separate list of np.arrays.
        :return: Either one list of np.arrays or a tuple of two such lists, as described in more detail above.
        """
        trajectories = []
        velocities = []

        with h5py.File(self.filename, self.mode) as h5f:
            dims = h5f.attrs['dimensions']
            for particle_id in h5f['particles']:
                if 'trajectory' in h5f[f'particles/{particle_id}'].keys():
                    trajectory = h5f[f'particles/{particle_id}/trajectory'][:]
                    trajectory_positions = trajectory[:, 1:dims + 1].transpose()
                    trajectory_velocities = trajectory[:, dims + 1:dims * 2 + 1].transpose()
                    trajectories.append(trajectory_positions)
                    if return_velocities:
                        velocities.append(trajectory_velocities)
                    # ax.plot(*trajectory_positions, color='grey')

        if return_velocities:
            return trajectories, velocities
        else:
            return trajectories


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End: