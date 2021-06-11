import numpy as np
import pytest

from cminject.fields.fluid_flow import StokesDragForceField


def test_stokes_field_inits_gas_params_from_file(example_field_file):
    field = StokesDragForceField(example_field_file)
    # the field file should be for nitrogen, expect appropriate values
    assert np.isclose(field.temperature, 293.15)
    assert np.isclose(field.dynamic_viscosity, 1.76e-5)
    assert np.isclose(field.m_gas, 4.7e-26)
    assert field.slip_correction_model == 'room_temp'


def test_stokes_field_can_have_gas_params_overridden(example_field_file):
    with pytest.warns(UserWarning):
        field = StokesDragForceField(example_field_file,
                                     temperature=273.15,
                                     dynamic_viscosity=4.2e-6,
                                     m_gas=12e-27)
    assert field.temperature == 273.15
    assert field.dynamic_viscosity == 4.2e-6
    assert field.m_gas == 12e-27


def test_stokes_field_warns_about_unknown_temps_and_uses_room_temp(example_field_file):
    # Expect a warning to trigger for 100.0K
    with pytest.warns(UserWarning) as record:
        field = StokesDragForceField(example_field_file,
                                     temperature=100.0)
        assert 'temperature' in str(record[0].message).lower()
        assert field.dynamic_viscosity == 1.76e-5

    # Expect that no warning happens for 293.15K
    with pytest.warns(None) as record:
        field = StokesDragForceField(example_field_file,
                                     temperature=293.15)
    assert len(record) == 0


def test_stokes_field_gives_reasonable_accelerations_for_example_file(example_field_file, spherical_particle_2d):
    field = StokesDragForceField(example_field_file)

    # choose some coordinate that's not on a grid point
    spherical_particle_2d.position[1] = 1.73e-5
    a = field.calculate_acceleration(spherical_particle_2d, time=0.0)
    assert np.isclose(a[0], 0.0)  # vr should be numerically zero
    assert np.isclose(a[1], 77.12521, atol=1e-5, rtol=1e-7)  # vz should be close to a known value


def test_stokes_field_returns_nan_vector_for_outside_particle(example_field_file, spherical_particle_2d):
    field = StokesDragForceField(example_field_file)

    # negative check: the particle position (0, 0) is inside and should *not* yield NaNs
    spherical_particle_2d.position = np.array([0.0, 0.0])
    a = field.calculate_acceleration(spherical_particle_2d, time=0.0)
    assert all(~np.isnan(a))

    # positive checks: positions we know to be outside the flow field should yield NaNs
    spherical_particle_2d.position = np.array([-1.0, 0.0])
    a = field.calculate_acceleration(spherical_particle_2d, time=0.0)
    assert all(np.isnan(a))

    spherical_particle_2d.position = np.array([0.0, 1.0])
    a = field.calculate_acceleration(spherical_particle_2d, time=0.0)
    assert all(np.isnan(a))

    # also verify consistency with the results of is_particle_inside
    randpos_x = np.random.normal(0, 1, 1000)
    randpos_y = np.random.normal(0, 1, 1000)
    for px, py in zip(randpos_x, randpos_y):
        spherical_particle_2d.position = np.array([px, py])
        if field.is_particle_inside(spherical_particle_2d.position, time=0.0):
            assert all(~np.isnan(a))
        else:
            assert all(np.isnan(a))
