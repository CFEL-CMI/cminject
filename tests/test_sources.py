import numpy as np

from cminject.particles import SphericalParticle
from cminject.sources import VariableDistributionSource
from cminject.utils.distributions import GaussianDistribution, UniformDistribution, NormOfDistributions, \
    DiracDeltaDistribution


def test_variable_dist_source_uses_given_distributions(dims2d):
    dist1 = GaussianDistribution(mu=0.0, sigma=1.0)
    dist2 = UniformDistribution(min=0.0, max=1.0)
    dist3 = NormOfDistributions(GaussianDistribution(mu=0.0, sigma=1.0),
                                GaussianDistribution(mu=0.0, sigma=1.0))
    dist4 = DiracDeltaDistribution(42.013)

    src = VariableDistributionSource(
        10000,
        position=[dist1, dist2],
        velocity=[dist2, dist1],
        radius=dist3,
        density=dist4,
        particle_class=SphericalParticle
    )
    particles = src.generate_particles()

    assert len(particles) == 10000  # sanity check
    assert np.isclose(np.std([p.position[0] for p in particles]), 1.0, atol=1e-2, rtol=1e-2)
    assert np.all([0 <= p.position[1] <= 1 for p in particles])
    assert np.isclose(np.std([p.velocity[1] for p in particles]), 1.0, atol=1e-2, rtol=1e-2)
    assert np.all([0 <= p.velocity[0] <= 1 for p in particles])
    assert np.all([p.radius > 0 for p in particles])
    assert np.all([p.rho == 42.013 for p in particles])
