import numpy as np
from cminject.base import Detector

from cminject.detectors import SimpleZDetector
from cminject.particles import SphericalParticle
from cminject.utils.global_config import GlobalConfig, ConfigKey


def test_detector_returning_true_updates_particle_detector_hits():
    class SillyDetector(Detector):
        def _has_particle_reached_detector(self, particle_identifier: int, position: np.array):
            return True
        def _hit_position(self, phase_space_position: np.array):
            return phase_space_position[:2] + 1.0  # ah well - it's just for testing
        def z_boundary(self):
            return (-np.inf, np.inf)

    det = SillyDetector('silly')
    p = SphericalParticle(0, 0, np.random.normal(0, 1, 2), np.random.normal(0, 1, 2),
                          radius=1e-6, rho=1050)
    assert len(p.detector_hits) == 0
    assert det.try_to_detect(p)
    assert len(p.detector_hits) == 1
    assert 'silly' in p.detector_hits
    assert np.isclose(p.detector_hits['silly'][0]['position'] - 1, p.position).all()
    assert (p.detector_hits['silly'][0]['velocity'] == p.velocity).all()


def test_simple_z_detector():
    """
    Tests the basic functionality of the SimpleZDetector class.
    """
    GlobalConfig().set(ConfigKey.NUMBER_OF_DIMENSIONS, 2)
    det = SimpleZDetector(0.0, identifier='test0')
    pos = np.array([1.0, 0.0])
    vel = np.array([1.0, 2.0])
    p = SphericalParticle(0, 0, pos, vel, radius=1e-6, rho=1050)

    # start exactly on detector -- should detect
    p.position[-1] = 0.0
    assert det.try_to_detect(p)
    assert 'test0' in p.detector_hits
    assert len(p.detector_hits['test0']) == 1
    hit = p.detector_hits['test0'][-1]
    # TODO why is there an additional dimension of length 1 here?
    assert hit['position'][0, 0] == 1.0
    assert hit['position'][0, 1] == 0.0

    # move to the right -- should not detect again
    p.position[-1] = 0.1
    assert not det.try_to_detect(p)
    assert len(p.detector_hits['test0']) == 1

    # move to the left crossing the detector -- should detect again
    p.position[-1] = -0.1
    assert det.try_to_detect(p)
    assert len(p.detector_hits['test0']) == 2

    # move to the right landing exactly on the detector -- should detect again
    p.position[-1] = 0.0
    assert det.try_to_detect(p)
    assert len(p.detector_hits['test0']) == 3

    # move further to the right -- should not detect again
    p.position[-1] = 0.1
    assert not det.try_to_detect(p)
    assert len(p.detector_hits['test0']) == 3
