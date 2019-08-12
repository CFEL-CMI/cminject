#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of CMInject
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

"""
Utility functions for handling commandline arguments.
"""

import argparse


def dist_description(x):
    """
    Defines a custom argparse type for a distribution description.
    :param x: The string value to parse.
    :return: A dictionary representing the parsed distribution description.
    :raises: argparse.ArgumentTypeError
    """
    shorthand_mapping = {
        'G': 'gaussian', 'RG': 'radial_gaussian',
        'L': 'linear', 'RL': 'radial_linear',
        'U': 'uniform', 'RU': 'radial_uniform'
    }
    values = x.split()
    kind = values[0]
    if kind in ['G', 'RG']:
        if len(values) != 3:
            raise argparse.ArgumentTypeError("Need a mu and sigma for a gaussian distribution description!")
        return {'kind': shorthand_mapping[kind],
                'mu': float(values[1]),
                'sigma': float(values[2])}
    elif kind in ['L', 'RL', 'U', 'RU']:
        if len(values) != 3:
            raise argparse.ArgumentTypeError("Need a min and max for a uniform distribution description!")
        return {'kind': shorthand_mapping[kind],
                'min': float(values[1]),
                'max': float(values[2])}
    else:
        raise argparse.ArgumentTypeError(f"Unknown distribution kind: {kind}!")


def natural_number(x):
    """
    Defines a custom argparse type for a natural number (excluding 0).
    :param x: The string value to parse.
    :return: An int value containing the natural number.
    :raises: argparse.ArgumentTypeError
    """
    x = int(x, base=10)
    if x < 1:
        raise argparse.ArgumentTypeError("A natural number was expected, minimum is 1!")
    return x


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
