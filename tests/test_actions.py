import numpy as np

from cminject.actions.brownian_motion import StokesBrownianMotionStep
from cminject.fields.fluid_flow import StokesDragForceField


def test_changes_particle_position_and_velocity(spherical_particle_2d, example_field_file,
                                                dims2d, dt1e_5):
    field = StokesDragForceField(example_field_file)
    bm_step = StokesBrownianMotionStep(field)

    orig_pos = np.array([0.0, 0.0])
    spherical_particle_2d.position = orig_pos
    bm_step(spherical_particle_2d, time=0.0)
    assert (orig_pos != spherical_particle_2d.position).any()


def test_doesnt_change_particle_position_and_velocity_when_outside(spherical_particle_2d, example_field_file,
                                                                   dims2d, dt1e_5):
    field = StokesDragForceField(example_field_file)
    bm_step = StokesBrownianMotionStep(field)

    orig_pos = np.array([0.0, -1.0])
    spherical_particle_2d.position = orig_pos
    bm_step(spherical_particle_2d, time=0.0)
    assert (orig_pos == spherical_particle_2d.position).all()
