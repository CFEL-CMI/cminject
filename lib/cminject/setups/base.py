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
from typing import List

from cminject.experiment import Experiment


class Setup(ABC):
    """
    A base class for classes that define experiment setups. Should be considered static, i.e. its method are all
    "@staticmethod"s, and no instances should (need to be) created.
    """
    @staticmethod
    @abstractmethod
    def parse_args(argarr: List[str]) -> argparse.Namespace:
        """
        Parse the arguments relevant to this specific setup. Will receive only the arguments that the main program
        (bin/cminject) did not recognise.
        :param argarr: The list of argument strings that were not recognised by the main program and so should be
            parsed by this setup class.
        :return: An argparse.Namespace, i.e. a result of calling parse.parse_args(argarr) for an appropriately
            constructed ArgumentParser instance.
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