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
import pathlib
import random
import warnings
from typing import Dict, Any, List, Tuple

import h5py
import numpy as np

from cminject.definitions.particles.base import Particle
from cminject.definitions.util import ParticleDetectorHit

from .base import ResultStorage


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
    dimensions used in the simulation. Further attributes can be stored on the root node by passing a metadata
    dictionary to the __init__ method.

    .. note::
        For efficient and predictable storage, it's assumed that all detector hits on each detector return
        exactly the same shape of data, and the same property description. This should be the case in any sensible
        implementation of a detector anyways. If you do find you need to deviate from this, store an array of the
        maximum necessary length for each hit and fill out values you can't have for a specific hit with `np.nan`.
    """
    def __init__(self, filename: str, mode: str = 'r', metadata: Dict[str, Any] = None):
        """
        Constructs a HDF5ResultStorage object.

        :param filename: The filename to store in/read from.
        :param mode: The mode to open the file with. 'r' by default to avoid overwriting data when using the convenience
            methods to read data from a written file, so if an instance should be used to actually store data,
            'w' or similar must be passed explicitly (refer to `h5py.File` docs for available modes).
        :param metadata: A dictionary of metadata to store on the output file for reference, like for instance
             all parameters involved in the experiment that generated the result.
        """
        self.filename = filename
        self.mode = mode
        self.metadata = metadata

        # Verify that the output file doesn't already exist, let's avoid overwriting existing data
        if self.mode != 'r' and os.path.isfile(filename):
            raise ValueError("Output file already exists! Please delete or move the existing file.")

        # Try to create the directories the file should be stored in
        pathlib.Path(os.path.dirname(self.filename)).mkdir(parents=True, exist_ok=True)

    def store_results(self, particles: List[Particle]) -> None:
        """
        Stores experiment results into an HDF5 file as described in the class docstring.

        :param particles: The list of particles as returned as the result from a finished Experiment run.
        :return: Nothing.
        """
        with h5py.File(self.filename, self.mode) as h5f:
            dimensions = len(particles[0].spatial_position)
            h5f.attrs['dimensions'] = dimensions

            # Store the metadata: Try storing each value directly, if that fails, try storing a string representation,
            # and if anything fails at that point, output a warning. We do NOT want to lose experiment results just
            # because metadata could not be stored.
            for k, v in self.metadata.items():
                try:
                    h5f.attrs[k] = v
                except TypeError:
                    try:
                        h5f.attrs[k] = str(v)
                    except Exception as e:
                        warnings.warn(f"Could not store metadata with key {k}, exception occurred: {e}")

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

    def get_dimensions(self) -> int:
        """
        Gets the number of spatial dimensions stored on the HDF5 file. Useful for properly separating out e.g. positions
        and velocities from the arrays returned by get_trajectories and get_positions.

        :return: The number of spatial dimensions.
        """
        with h5py.File(self.filename, self.mode) as h5f:
            return int(h5f.attrs['dimensions'])

    def get_trajectories(self, n_samples: int = None, positive_r: bool = False) -> List[np.array]:
        """
        Gets all stored trajectories as a list of np.arrays, each array being (d, n_i)-shaped, where d is the number
        of measured quantities of the stored results, and n_i is the number of data points along the trajectory
        of the i-th particle.

        :param n_samples: The number of samples to draw (without replacement) from the whole list of trajectories.
            If None, all trajectories are returned. Note that if there are less than n_samples trajectories in the file,
            this will not fail, but all trajectories in the file will be returned.
        :return: A list of np.arrays as described above.
        """
        trajectories = []

        with h5py.File(self.filename, self.mode) as h5f:
            keys = list(h5f['particles'].keys())
            if n_samples is not None:
                n_trajs = len(keys)
                ids = random.sample(keys, min(n_samples, n_trajs))
            else:
                ids = keys

            for particle_id in ids:
                if 'trajectory' in h5f[f'particles/{particle_id}'].keys():
                    trajectory = h5f[f'particles/{particle_id}/trajectory'][:].transpose()
                    if positive_r:
                        trajectory[1] = np.abs(trajectory[1])
                    trajectories.append(trajectory)

        return trajectories

    def get_trajectories_iterator(self, n_samples: int = None, abs_indices: List[int] = None):
        """
        Gets all stored trajectories as a list of np.arrays, each array being (d, n_i)-shaped, where d is the number
        of measured quantities of the stored results, and n_i is the number of data points along the trajectory
        of the i-th particle.

        :param n_samples: The number of samples to draw (without replacement) from the whole list of trajectories.
            If None, all trajectories are returned. Note that if there are less than n_samples trajectories in the file,
            this will not fail, but all trajectories in the file will be returned.
        :return: A list of np.arrays as described above.
        """
        if abs_indices is None:
            abs_indices = []

        with h5py.File(self.filename, self.mode) as h5f:
            keys = list(h5f['particles'].keys())
            if n_samples is not None:
                n_trajs = len(keys)
                ids = random.sample(keys, min(n_samples, n_trajs))
            else:
                ids = keys

            for particle_id in ids:
                if 'trajectory' in h5f[f'particles/{particle_id}'].keys():
                    trajectory = h5f[f'particles/{particle_id}/trajectory'][:].transpose()
                    for abs_index in abs_indices:
                        trajectory[abs_index] = np.abs(trajectory[abs_index])
                    yield trajectory

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

    def get_detectors(self) -> List[Tuple[Any, np.array]]:
        """
        Gets all stored detectors with their hits.

        :return: A list of 2-tuples, the first component being each detector's identifier and the second
            being a (d, m_i)-shaped np.array containing all hits.
        """
        detectors = []
        with h5py.File(self.filename, self.mode) as h5f:
            if 'detector_hits' in h5f:
                for detector_id in sorted(h5f['detector_hits'], key=int):
                    hits = h5f['detector_hits'][detector_id][:].transpose()
                    detectors.append((detector_id, hits))
        return detectors

    def get_reconstructed_detectors(self) -> List[Tuple[float, np.array]]:
        """
        Gets all reconstructed detectors as stored by the `cminject_reconstruct-detectors` tool.

        :return: A list of 2-tuples, the first component being each detector's Z position and the second
            being a (d, k_i)-shaped np.array containing all hits, where k_i is the number of reconstructed quantities
            (the number of --xi parameters passed to `cminject_reconstruct-detectors`).
        """
        with h5py.File(self.filename, self.mode) as h5f:
            if 'reconstructed_detectors' in h5f:
                dataset = h5f['reconstructed_detectors']
                return [(float(z), dataset[z][:]) for z in dataset.keys()]
            else:
                return []

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
