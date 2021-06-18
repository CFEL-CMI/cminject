"""
These tests are essentially full-software integration tests. They function by
calling the example simulation scripts in the `examples/` folder, then making
assertions about the results of these simulations.
"""

import os
import subprocess
import pathlib
import tempfile
import numpy as np

import h5py
import pytest

from cminject.result_storages import HDF5ResultStorage

@pytest.fixture
def examples_path():
    return pathlib.Path(__file__).parent.absolute() / '..' / '..' / 'examples'


def test_example_simulation(examples_path):
    """
    Verifies the overall functionality of the simple example setup. It does so
    by running a simulation with 100 particles, letting it write to a temporary
    output HDF5 file, and making some assertions about the basic structure of
    this file, by opening it via HDF5ResultStorage:

    - the keys 'particles' and 'detectors' should be present
    - there should be 100 trajectories stored
    - there should be 3 detectors stored
    - there should be a detector at z=0.0, with at least one particle hit
    """
    nparticles = 100

    os.chdir(examples_path)
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, 'example_output.h5')
        output = subprocess.check_output(f"""
        cminject -s cminject.setups.example.ExampleSetup \
            -f example_field.h5 \
            -n {nparticles} \
            -o {outfile} \
            -T
        """, shell=True)

        with h5py.File(outfile, 'r') as h5f:
            assert 'particles' in h5f.keys()
            assert 'detectors' in h5f.keys()

        with HDF5ResultStorage(outfile, 'r') as storage:
            assert storage.get_dimensions() == 2

            trajectories = storage.get_trajectories()[:]
            assert trajectories.shape == (nparticles,)
            assert any(len(trajectory) > 0 for trajectory in trajectories)

            detectors = storage.get_detectors()
            assert len(list(detectors.keys())) == 3
            assert 'SimpleZ@0.0' in detectors.keys()
            assert len(detectors['SimpleZ@0.0'][:]) > 0


def test_als_simulation(examples_path):
    """
    Verifies that results of an example ALS simulation are reasonable output.
    This means that particle hits are present at certain detectors, and that
    they have a distribution with mean and standard deviation within a
    reasonable range.

    These tests are somewhat lenient in these statistical checks, to avoid
    erroring due to low statistics or chance. They are only intended to verify
    that the software behaves as expected, not to provide a guarantee that an
    experiment is well-modeled or anything in that direction.
    """
    nparticles = 500

    os.chdir(examples_path)
    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, 'als_output.h5')
        output = subprocess.check_output(f"""
        cminject -f "als_field.h5"\
            -rho 19320\
            -D 2\
            -p G[0,0.002] -0.128\
            -v G[0,1] G[1,0.1]\
            -r G[1.35e-8,1.91e-9]\
            -t 0 1\
            -ts 1e-5\
            -d 0.001 0.00141 0.00181 0.00221 0.00261\
               0.00301 0.00341 0.00381 0.00421 0.00461\
            -n {nparticles}\
            -o {outfile}\
            -B \
            -T
        """, shell=True)

        with h5py.File(outfile, 'r') as h5f:
            assert 'particles' in h5f.keys()
            assert 'detectors' in h5f.keys()

        with HDF5ResultStorage(outfile, 'r') as storage:
            assert storage.get_dimensions() == 2

            trajectories = storage.get_trajectories()[:]
            assert trajectories.shape == (nparticles,)
            assert any(len(trajectory) > 0 for trajectory in trajectories)

            # at least half of the particles must have crossed the exit
            nparticles_after_exit = sum(any(trajectory['position'][:, 1] > 0.0) for trajectory in trajectories)
            assert nparticles_after_exit > nparticles//2

            detectors = storage.get_detectors()
            assert len(list(detectors.keys())) == 10
            assert 'SimpleZ@0.00301' in detectors.keys()
            hits = detectors['SimpleZ@0.00301'][:]
            assert len(hits) > nparticles//2
            r = hits['position'][:, 0]
            mean_r = np.mean(r)
            std_r = np.std(r)
            focus_size = 4 * np.percentile(r, 70)
            assert 25e-6 <= focus_size <= 50e-6
            assert np.abs(mean_r) < 5e-6
            assert 0 <= std_r <= 3e-5
