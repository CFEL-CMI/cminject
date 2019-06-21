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

    NOTE: For efficient and predictable storage, it's assumed that all detector hits on each detector return
    exactly the same shape of data, and the same property description. This should be the case in any sensible
    implementation of a detector anyways. If you do find you need to deviate from this, store an array of the maximum
    necessary length for each hit and fill out values you can't have for a specific hit with `np.nan`.
    """
    def __init__(self, filename: str):
        self.filename = filename

    def store_results(self, particles: List[Particle]) -> None:
        with h5py.File(self.filename) as h5f:
            detector_hits = {}

            for particle in particles:
                if particle.trajectory:
                    h5f[f'trajectories/{particle.identifier}'] = np.array(particle.trajectory)
                    h5f[f'trajectories/{particle.identifier}'].attrs['description'] = \
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
