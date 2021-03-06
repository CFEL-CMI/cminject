#!/usr/bin/env python3
# We need argparse for parsing arguments from the command line
import argparse

# Import the most basic definitions for defining a Setup which returns an Experiment
# and has arguments parsed by a SetupArgumentParser
from cminject.experiment import Experiment
from cminject.base import Setup
from cminject.utils.args import SetupArgumentParser, distribution_description

# Import some concrete class implementations to be able to define a simple setup:
from cminject.utils.distributions import constant, GaussianDistribution  # Required types of distributions
from cminject.sources import VariableDistributionSource  # Generates particles from property distributions
from cminject.detectors import SimpleZDetector  # A detector which is positioned at some Z position
# A device for fluid flow based on an HDF5 file, together with FlowType, an enumeration of the possible models to use
from cminject.devices.fluid_flow import FluidFlowDevice, FlowType


class ExampleSetup(Setup):
    """
    A simple setup which will simulate a single flow field with particles starting from a fixed position and velocity
    distribution and radius. The only parameters available are the flow field filename and the density (density) of the
    particle material. The density is there just to show how additional arguments can be added.
    """
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        dt = 1e-5
        experiment = Experiment(number_of_dimensions=2, time_step=dt, time_interval=(0, 1))
        experiment.add_source(VariableDistributionSource(
            number_of_particles=main_args.nof_particles,
            position=args.pos,
            velocity=[GaussianDistribution(mu=0.0, sigma=1e-3), GaussianDistribution(mu=2, sigma=0.5)],
            radius=constant(5e-9), density=constant(args.rho)
        ))
        experiment.add_device(FluidFlowDevice(filename=args.filename, flow_type=FlowType.STOKES, brownian_motion=True))
        for z in [0.0, 0.01, 0.02]:
            experiment.add_detector(SimpleZDetector(z))

        return experiment

    @staticmethod
    def get_parser() -> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-f', '--filename', help='The filename of the flow field (HDF5).', type=str)
        parser.add_argument('--rho', help='The density of the particle material [kg/m^3].', type=float, default=1050.0)
        parser.add_argument('--pos', help='The position distributions in x/z space [m]', type=distribution_description,
                            nargs=2, default=[GaussianDistribution(0, 1e-3), constant(0.0)])
        return parser
