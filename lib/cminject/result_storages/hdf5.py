#!/usr/bin/env python3
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

"""
Implementations of :class:`cminject.base.ResultStorage` that store results as HDF5 files.
"""

import os
import warnings
from typing import Dict, Any, List, Union

import h5py
import numpy as np

from cminject.base import Particle, ResultStorage, ParticleDetectorHit


class HDF5ResultStorage(ResultStorage):
    """
    Stores Experiment output in an HDF5 file.

    Let n be the number of stored particles, d_p be the number of dimensions in the particles' `position` property,
    d_t be the number of dimensions in a point along a trajectory, and d_d the number of dimensions stored on each
    `ParticleDetectorHit`. All d_* variables depend on the concrete implementations used for a simulation.

    The following groups and datasets will be stored:

    # FIXME out of date
    - particles:

      - initial_positions: The initial phase space positions of all particles, as one (n, d_p) array.
      - final_position: The final phase space positions of all particles, as one (n, d_p) array.
      - trajectories: A jagged array of (n, d_t, X) entries, where X can be variable (jagged)..
        This dataset is not stored if no nonempty trajectories were tracked.

    - detectors/{id}: All hits on the detector with identifier <id>, as one (h_i, d_d) array. h_i is the number of hits
      on the i-th detector. Each entry is one hit, one particle can incur multiple hits, and what is stored at each
      hit is the concatenation of (position, properties) of the detected particle.

    .. note::
      Particle trajectories are stored as one `awkward.JaggedArray` (see https://github.com/scikit-hep/awkward-array),
      meaning that when reading from such a result file, you either need to use this class to read out the results,
      or handle this stored format yourself.

    Further attributes can be stored on the root node by passing a dictionary as the "metadata" argument to the
    constructor.

    .. note::
      For efficient and predictable storage, it's assumed that all detector hits on each detector return
      exactly the same shape of data, and the same property description. This should probably be the case in any
      a typical implementation of a detector.
    """
    def __init__(self, filename: str, mode: str = 'r', force_writable: bool = False, metadata: Dict[str, Any] = None):
        self.filename = filename
        self.mode = mode
        self.metadata = metadata
        self._handle: Union[h5py.File, None] = None  # will be set in __enter__

        if self.mode != 'r' and not force_writable:
            # Verify that the output file doesn't already exist, let's avoid overwriting existing data
            if os.path.isfile(filename):
                raise ValueError(f"Output file {filename} already exists! Please delete or move the existing file.")

    @property
    def file_handle(self):
        """
        Returns the file handle (a h5py.File instance) that this result storage currently has opened.

        :return: The file handle.
        :raises: ``RuntimeError``, if the file handle is not set (e.g., due to the result storage being accessed outside
          a with-block)
        """
        return self._require_handle()

    def _require_handle(self):
        if getattr(self, '_handle', None) is None:
            raise RuntimeError(
                f"File handle is not set! {self.__class__.__name__} must be used as a with-block: \n"
                f"    with {self.__class__.__name__}(...) as storage:\n"
                f"        storage.some_method()\n"
                f"        ..."
            )
        return self._handle

    def __enter__(self):
        self._handle = h5py.File(self.filename, self.mode)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._handle.close()
        self._handle = None

    def store_results(self, particles: List[Particle]) -> None:
        if not particles:
            return

        f = self._require_handle()
        dims = particles[0].number_of_dimensions
        f.attrs['dimensions'] = dims
        # Store the metadata: Try storing each value directly, if that fails, try storing a string representation,
        # and if anything fails at that point, output a warning. We do NOT want to lose experiment results just
        # because metadata could not be stored.
        for k, v in self.metadata.items():
            try:
                f.attrs[k] = v
            except TypeError:
                try:
                    f.attrs[k] = str(v)
                except Exception as e:
                    warnings.warn(f"Could not store metadata with key {k}, exception occurred: {e}")

        # Get the appropriate dtypes, assuming all particles are instances of the same class
        const_dtype = particles[0].as_dtype('constant')
        tracked_dtype = particles[0].as_dtype('tracked')
        all_dtype = particles[0].as_dtype('all')

        # Generate arrays that can be efficiently stored and accessed
        n = len(particles)
        pg = f.create_group('particles')
        props = pg.create_dataset('properties', (n,), dtype=const_dtype)

        tg = pg.create_group('tracked')
        tracked_initial = tg.create_dataset('initial', (n,), dtype=tracked_dtype)
        tracked_final = tg.create_dataset('final', (n,), dtype=tracked_dtype)

        trajectories = []

        # Aggregate data into arrays prepared above
        detector_hits: Dict[Any, List[ParticleDetectorHit]] = {}
        for i, particle in enumerate(particles):
            props[i] = particle.as_array('constant')
            tracked_initial[i] = particle.initial_tracked_properties
            tracked_final[i] = particle.as_array('tracked')

            trajectory = particle.trajectory or []
            trajectories.append(trajectory)

            if particle.detector_hits:
                for detector_id, hits in particle.detector_hits.items():
                    if detector_id not in detector_hits:
                        detector_hits[detector_id] = []
                    detector_hits[detector_id] += hits

        # Store trajectories
        if np.any(len(traj) for traj in trajectories):  # if any of the particles had a non-empty trajectory
            vlen_dtype = h5py.vlen_dtype(tracked_dtype)
            ds = tg.create_dataset('trajectories', (n,), dtype=vlen_dtype)
            for i in range(n):
                ds[i] = np.array(trajectories[i], dtype=tracked_dtype).ravel()

        # Store detectors
        dg = f.create_group('detectors')
        for detector_id, hits in detector_hits.items():
            if hits:
                dg.create_dataset(str(detector_id), (len(hits),), dtype=all_dtype,
                                  data=np.array([hit.hit_state for hit in hits]))

    def get_dimensions(self) -> int:
        f = self._require_handle()
        return int(f.attrs['dimensions'])

    def get_properties(self) -> np.ndarray:
        f = self._require_handle()
        return f.get('particles/properties')

    def get_trajectories(self) -> np.ndarray:
        f = self._require_handle()
        return f.get('particles/tracked/trajectories')

    def get_tracked_initial(self) -> np.ndarray:
        f = self._require_handle()
        return f.get('particles/tracked/initial')

    def get_tracked_final(self) -> np.ndarray:
        f = self._require_handle()
        return f.get('particles/tracked/final')

    def get_detectors(self) -> Dict[str, np.ndarray]:
        f = self._require_handle()
        detectors = {}
        if 'detectors' in f:
            for detector_id in f['detectors']:
                hits = f['detectors'][detector_id][:].T
                detectors[detector_id] = hits
        return detectors or None

    # TODO store description as metadata!!
    # TODO maybe add function to return a pandas.DataFrame (?)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
