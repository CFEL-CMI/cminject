import numpy as np
import argparse

from cminject.definitions.devices.stark_deflector import StarkDeflector
from cminject.definitions.particles import Molecule

from cminject.definitions.setups import Setup
from cminject.definitions.sources import MolDistributionSource
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
            z_minmax=args.z_minmax, field_limit=field_limit)
        ]

        detectors = [SimpleZDetector(identifier=i, z_position=pos) for i, pos in enumerate(args.detectors)]
        sources = [MolDistributionSource(main_args.nof_particles, position=args.position, velocity=args.velocity, temperature=temperature,
                                         jmax = Jmax, mass=args.mass, energy_filename=energy_filename)]


        return Experiment(devices=devices, detectors=detectors, sources=sources, property_updaters=property_updaters,
                          time_interval=(t_start, t_end), time_step=dt, seed=main_args.seed,
                          number_of_dimensions=args.dimensions)

    @staticmethod
    def get_parser()-> SetupArgumentParser:
        parser = SetupArgumentParser()
        parser.add_argument('-f', '--field-Gradient-filename', help='the norm and gradient of the Stark field (hdf5 format)', type=str, required=True)

        parser.add_argument('-fl', '--field-limit', help='Maximum value of the norm of the electric field to be considered as a point on an electrod',
                            type=int, required=True)  # should be two or three? should I omit it?

        parser.add_argument('-enr', '--energy-filename', help='Energies of the molecule as a function of the norm of the field (hdf5 format)',
                            type=int, required=True)  # should be two or three? should I omit it?

        parser.add_argument('-D', '--dimensions', help='# of spatial dimensions',
                            type=int, required=True)  # should be two or three? should I omit it?

        parser.add_argument('-p', '--position',
                            help='Distribution description for the position.',
                            nargs='*', type=dist_description, required=True)

        parser.add_argument('-v', '--velocity',
                            help='Distribution description for the velocity.',
                            nargs='*', type=dist_description, required=True)

        parser.add_argument('-T', '--temperature',
                            help='Tempreture of the molecular beam in Kelvin.',
                            type=int, required=True)
        """
        where and how to include this piece of info?
        """
        parser.add_argument('-j', '--Jmax',
                            help='Up to what Jmax shall the simulation be carried?',
                            type=int, required=True)

        parser.add_argument('-z', '--z-minmax',
                            help='z_start and z_final of the z-dimension',
                            type=Tuple[float, float], required=True)

        parser.add_argument('-d', '--detectors', help='The Z positions of the detectors', nargs='+',
                            type=float, required=True)

        return parser
