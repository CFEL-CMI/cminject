import numpy as np

from cminject.boundaries.field_based import GridFieldBasedBoundary
from cminject.boundaries.simple import SimpleZBoundary, InfiniteBoundary,\
    CuboidBoundary
from cminject.fields.fluid_flow import StokesDragForceField


def test_simple_z_boundary(dims2d):
    boundary = SimpleZBoundary(z_minmax=(-0.1, 0.9))
    positions_and_results = [
        ([0.0, 0.0], True),
        ([0.0, -0.1], True),
        ([0.0, 0.9], True),
        ([0.0, -0.100001], False),
        ([0.0, 0.900001], False),
        ([1e-12, 0.0], True),
        ([1e80, -0.1], True),
        ([np.nan, 0.9], True),
        ([-11891111111111, -0.100001], False),
        ([np.nan, 0.900001], False)
    ]

    for (position, expected_result) in positions_and_results:
        result = boundary.is_particle_inside(np.array(position, dtype=np.float64), time=0.0)
        errorstr = f"{position} should be {'inside' if expected_result else 'outside'} boundary but isn't"
        assert result == expected_result, errorstr


def test_cuboid_boundary(dims2d):
    boundary = CuboidBoundary(intervals=[(-0.5, 0.1), (-0.1, 0.9)])
    positions_and_results = [
        ([0.0, 0.0], True),
        ([0.0, -0.1], True),
        ([0.0, 0.9], True),
        ([-0.6, 0.0], False),
        ([-0.5000001, -0.1], False),
        ([0.1232, 0.9], False),
        ([0.0, -0.100001], False),
        ([0.0, 0.900001], False),
        ([1e-12, 0.0], True),
        ([1e80, -0.1], False),
        ([np.nan, 0.9], False),
        ([np.nan, np.nan], False),
        ([-11891111111111, -0.100001], False),
        ([np.nan, 0.900001], False)
    ]

    for (position, expected_result) in positions_and_results:
        result = boundary.is_particle_inside(np.array(position, dtype=np.float64), time=0.0)
        errorstr = f"{position} should be {'inside' if expected_result else 'outside'} boundary but isn't"
        assert result == expected_result, errorstr


def test_infinite_boundary(dims2d):
    boundary = InfiniteBoundary()

    positions_and_results = [
        ([0.0, 0.0], True),
        (np.random.normal(0, 1e130, 2), True),
        ([np.nan, np.nan], True),
        (np.random.normal(0, 1e300, 2), True),
        ([-np.inf, np.inf], True)
    ]

    for (position, expected_result) in positions_and_results:
        result = boundary.is_particle_inside(np.array(position, dtype=np.float64), time=0.0)
        errorstr = f"{position} should be {'inside' if expected_result else 'outside'} boundary but isn't"
        assert result == expected_result, errorstr


def test_grid_field_based_boundary_outside_interpolates_to_nan(dims2d, example_field_file):
    field = StokesDragForceField(example_field_file)
    positions = np.random.normal(0, 1, (10000, 2))
    boundary = GridFieldBasedBoundary(field)

    for position in positions:
        field_says_inside = field.is_particle_inside(position, time=0.0)
        boundary_says_inside = boundary.is_particle_inside(position, time=0.0)
        assert field_says_inside == boundary_says_inside
