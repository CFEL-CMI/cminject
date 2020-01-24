import numpy as np
import argparse

from cminject.definitions.devices.Stark_deflector import StarkDeflector
from cminject.definitions.particles import Molecule

from cminject.definitions.setups import Setup
from cminject.definitions.sources import VariableDistributionSource
from cminject.experiment import Experiment
from cminject.utils.args import dist_description, SetupArgumentParser


class StarkExp(Setup):
    """
    A setup for deflection molecular beams in a Stark deflector.
    """
    @staticmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        t_start, t_end = main_args.time_interval
        dt = main_args.time_step

        devices = [StarkDeflector(
            filename=args.field_Gadient_filename,
            z_minmax=args.z_minmax)
        )]


        #detectors = [SimpleZDetector(identifier=i, z_position=pos) for i, pos in enumerate(args.detectors)]
        #sources = [VariableDistributionSource(main_args.nof_particles, position=args.position, velocity=args.velocity,
                                              radius=args.radius, rho=args.density, subclass=SphericalParticle)]

        return Experiment(devices=devices, detectors=detectors, sources=sources, property_updaters=property_updaters,
                          time_interval=(t_start, t_end), time_step=dt, seed=main_args.seed,
                          number_of_dimensions=args.dimensions)

    @staticmethod
    def get_parser()-> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-f', '--field-Gradient-filename', help='the norm and gradient of the Stark field (hdf5 format)', type=str, required=True)

        parser.add_argument('-D', '--dimensions', help='# of spatial dimensions',
                            type=int, required=True)

        parser.add_argument('-p', '--position',
                            help='Distribution description for the position.',
                            nargs='*', type=dist_description, required=True)

        parser.add_argument('-v', '--velocity',
                            help='Distribution description for the velocity.',
                            nargs='*', type=dist_description, required=True)
        parser.add_argument('-T', '--temperature',
                            help='Tempreture of the molecular beam in Kelvin.',
                            type=int, required=True)

        parser.add_argument('-j', '--Jmax',
                            help='Up to what Jmax shall the simulation be carried?',
                            type=int, required=True)

        parser.add_argument('-z', '--z-minmax',
                            help='z_start and z_final of the z-dimension',
                            type=Tuple[float, float], required=True)


        return parser
