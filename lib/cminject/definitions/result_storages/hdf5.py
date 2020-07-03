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

import os
import pathlib
import random
import warnings
from typing import Dict, Any, List, Tuple, Union

import h5py
import numpy as np
import awkward

from cminject.definitions.particles.base import Particle
from cminject.definitions.util import ParticleDetectorHit

from .base import ResultStorage


class HDF5ResultStorage(ResultStorage):
    """
    Stores Experiment output in an HDF5 file.

    Let n be the number of stored particles, d_p be the number of dimensions in the particles' `position` property,
    d_t be the number of dimensions in a point along a trajectory, and d_d the number of dimensions stored on each
    `ParticleDetectorHit`. All d_* variables depend on the concrete implementations used for a simulation.

    The following groups and datasets will be stored:

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
    def __init__(self, filename: str, mode: str = 'r', metadata: Dict[str, Any] = None):
        self.filename = filename
        self.mode = mode
        self.metadata = metadata
        self._handle: Union[h5py.File, None] = None  # will be set in __enter__

        if self.mode != 'r':
            # Verify that the output file doesn't already exist, let's avoid overwriting existing data
            if os.path.isfile(filename):
                raise ValueError(f"Output file {filename} already exists! Please delete or move the existing file.")

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
        f = self._require_handle()
        dims = len(particles[0].spatial_position)
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

        # Generate arrays that can be efficiently stored and accessed
        sdtype = h5py.string_dtype()  # special variable-length string dtype, for identifiers
        d = len(particles[0].initial_position)
        n = len(particles)
        pg = f.create_group('particles')
        ids = pg.create_dataset('identifiers', (n,), dtype=sdtype)
        initial_positions = pg.create_dataset('initial_positions', (n, d), dtype=np.float64)
        final_positions = pg.create_dataset('final_positions', (n, d), dtype=np.float64)
        # For trajectories, use a JaggedFillable to iteratively append to a jagged array that will be stored
        trajectories_fillable = awkward.generate.JaggedFillable(awkward.util.awkwardlib(None))

        # Aggregate data into arrays prepared above
        detector_hits: Dict[Any, List[ParticleDetectorHit]] = {}
        for i, particle in enumerate(particles):
            ids[i] = str(particle.identifier)
            initial_positions[i] = particle.initial_position
            final_positions[i] = particle.position
            trajectory = particle.trajectory or []
            trajectories_fillable.append(trajectory, awkward.generate.typeof(trajectory))

            if particle.detector_hits:
                for detector_id, hits in particle.detector_hits.items():
                    if detector_id not in detector_hits:
                        detector_hits[detector_id] = []
                    detector_hits[detector_id] += hits

        # Store trajectories
        trajectories = trajectories_fillable.finalize()
        if np.any(trajectories.starts != trajectories.stops):  # if any of the particles had a non-empty trajectory
            f_awkward = awkward.hdf5(f)
            f_awkward['particles/trajectories'] = trajectories
            f['particles/trajectories'].attrs['description'] = particles[0].position_description + ['t']

        # Store detectors
        dg = f.create_group('detectors')
        for detector_id, hits in detector_hits.items():
            if hits:
                dg[str(detector_id)] = np.array([hit.full_properties for hit in hits])
                dg[str(detector_id)].attrs['description'] = hits[0].full_properties_description

    def get_dimensions(self) -> int:
        f = self._require_handle()
        return int(f.attrs['dimensions'])

    def get_identifiers(self) -> np.ndarray:
        f = self._require_handle()
        return f.get('particles/identifiers')

    def get_trajectories(self) -> awkward.JaggedArray:
        f = self._require_handle()
        f_awkward = awkward.hdf5(f)
        return f_awkward.get('particles/trajectories')

    def get_initial_positions(self) -> np.ndarray:
        f = self._require_handle()
        return f.get('particles/initial_positions')

    def get_final_positions(self) -> np.ndarray:
        f = self._require_handle()
        return f.get('particles/final_positions')

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
