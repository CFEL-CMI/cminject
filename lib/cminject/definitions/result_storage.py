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

from typing import List

import h5py
import numpy as np

from .base import ResultStorage, Particle


class HDF5ResultStorage(ResultStorage):
    """
    Stores the results of an Experiment in an HDF5 file, with the following stored:

    - trajectories:

      - {i_0, i_1, ..., i_n}: The trajectory of the particle with identifier i.

    - detector_hits:

      - {j_0, j_1, ..., j_m}: All hits on the detector j. What's actually stored here is an
        array of the `full_properties` field of each `ParticleDetectorHit`.

    Each entry as mentioned above also has an attached `description` attribute which matches
    the stored array in length (if the class implementations have been done according to documentation),
    and describes/identifies the value at each position in the array.

    .. note::
        For efficient and predictable storage, it's assumed that all detector hits on each detector return
        exactly the same shape of data, and the same property description. This should be the case in any sensible
        implementation of a detector anyways. If you do find you need to deviate from this, store an array of the
        maximum necessary length for each hit and fill out values you can't have for a specific hit with `np.nan`.
    """
    def __init__(self, filename: str):
        self.filename = filename

    def store_results(self, particles: List[Particle]) -> None:
        with h5py.File(self.filename, 'w') as h5f:
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


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""