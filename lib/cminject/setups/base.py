#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of CMInject.
# It is a reference implementation of a subset of possible experimental setups, together with command line argument
# parsing. Its purpose is to allow those simulations to be run without additional code needing to be written, and to
# showcase how the setup is created so people can extend it or write more complex setups.
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this program for scientific work, you should correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

import argparse
from abc import ABC, abstractmethod

from cminject.experiment import Experiment
from cminject.utils.args import SetupArgumentParser


class Setup(ABC):
    """
    A base class for classes that define experiment setups. Should be considered static, i.e. its method are all
    "@staticmethod"s, and no instances should (need to be) created.
    """
    @staticmethod
    @abstractmethod
    def get_parser() -> SetupArgumentParser:
        """
        Returns a parser for the arguments relevant to this specific setup.
        This parser will receive only the arguments that the main program (bin/cminject) did not recognise.
        :return: An argparse.ArgumentParser (or subclass) instance that can parse all the args relevant to this setup.
        """
        pass

    @staticmethod
    def validate_args(args: argparse.Namespace):
        """
        Validates the arguments. Useful for validation that needs to check multiple args at once, and not just one
        specific arg. Overriding this method is optional and only needs to be done if such validation is desired.

        Should raise argparse.ArgumentError if validation failed, and have no effect if validation succeeded.
        :param args: The argparse.Namespace object that the parser constructed by get_parser() returned after being
            given all arguments that the main program did not recognise.
        :raises: argparse.ArgumentError
        """
        pass

    @staticmethod
    @abstractmethod
    def construct_experiment(main_args: argparse.Namespace, args: argparse.Namespace) -> Experiment:
        """
        Constructs an Experiment instance from the arguments to the main program and this specific setup.
        :param main_args: The argparse.Namespace parsed by the main program. Contains args that are common
            to all setups, like the number of particles and timespan.
        :param args: The argparse.Namespace parsed by this class's parse_args method.
        :return: An Experiment instance that's ready to be ran.
        """
        pass


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""