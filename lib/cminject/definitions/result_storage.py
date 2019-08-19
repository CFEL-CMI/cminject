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
from typing import List, Union, Tuple, Dict, Any

import h5py
import numpy as np

from .base import ResultStorage, Particle, ParticleDetectorHit


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
            detector_hits: Dict[Any, List[ParticleDetectorHit]] = {}

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


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End: