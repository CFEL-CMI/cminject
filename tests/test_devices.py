from cminject.actions.brownian_motion import StokesBrownianMotionStep, MolecularFlowBrownianMotionStep
from cminject.devices.fluid_flow import FluidFlowDevice, FlowType


def test_fluid_flow_field_device_adds_correct_brownian_motion_action_when_instructed(example_field_file):
    # Correct action for Stokes flow
    dev = FluidFlowDevice(example_field_file, flow_type=FlowType.STOKES,
                          brownian_motion=True)
    assert len(dev.actions) > 0
    assert any(
        isinstance(action, StokesBrownianMotionStep)
        for action in dev.actions
    ), "Field should contain a StokesBrownianMotionStep"

    # Correct action for molecular flow
    dev = FluidFlowDevice(example_field_file, flow_type=FlowType.MOLECULAR_FLOW,
                          brownian_motion=True)
    assert len(dev.actions) > 0
    assert any(
        isinstance(action, MolecularFlowBrownianMotionStep)
        for action in dev.actions
    ), "Field should contain a MolecularFlowBrownianMotionStep"

    # No action should be added when brownian_motion=False
    dev1 = FluidFlowDevice(example_field_file, flow_type=FlowType.STOKES,
                           brownian_motion=False)
    dev2 = FluidFlowDevice(example_field_file, flow_type=FlowType.MOLECULAR_FLOW,
                           brownian_motion=False)
    for dev in (dev1, dev2):
        assert not any(
            isinstance(
                action,
                (StokesBrownianMotionStep, MolecularFlowBrownianMotionStep)
            )
            for action in dev.actions
        ), "Field should not contain any kind of BrownianMotionStep"
