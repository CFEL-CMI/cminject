#!/usr/bin/env python3
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
import re
import warnings
from typing import Callable, Any

import numpy as np


class SetupHelpFormatter(argparse.MetavarTypeHelpFormatter):
    """
    A help formatter for argparse.ArgumentParsers that are returned by a Setup subclass's get_parser() method.
    Omits the "usage: <progname> arg1 [arg2] ..." line at the start.
    """
    def _format_usage(self, usage, actions, groups, prefix):
        return ""


class SetupArgumentParser(argparse.ArgumentParser):
    """
    An argparse.ArgumentParser subclass specifically suited to be returned by a Setup subclass's get_parser() method.
    Doesn't print that it has its own -h option (since that's never exposed), and uses the SetupHelpFormatter to
    avoid printing the usage.
    """
    def __init__(self, *args, **kwargs):
        if 'add_help' not in kwargs:
            kwargs['add_help'] = False
        if 'formatter_class' not in kwargs:
            kwargs['formatter_class'] = SetupHelpFormatter
        super().__init__(*args, **kwargs)


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
    if len(values) == 1:
        try:
            return float(values[0])
        except ValueError:
            raise argparse.ArgumentTypeError("When passing only one value, it must be parsable as a float!")

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


def parse_single_dimension_description(s: str) -> Callable[[np.array], Any]:
    """
    Parses dimension descriptions that can be of the following format:

    - x, y, z, vx, vy, vz refer to the first (x), second (y) and last (z) component of the position and velocity
      components
    - p<i> refers to the i-th entry of the position component
    - v<i> refers to the i-th entry of the velocity component
    - <some_name><i> refers to the i-th entry of the component called some_name
    - <some_name> refers to the component called some_name as a whole

    :param s: The string to parse
    :return: A Callable that will return the appropriate component from a given np.array, when given an array
        with a structured dtype that contains the given component. This Callable will throw an error when called if it
        can not execute the access the user wants on the given np.array.
    """
    s = s.strip()
    static = {
        'x': lambda x: x['position'][:, 0],
        'y': lambda x: x['position'][:, 1],
        'z': lambda x: x['position'][:, -1],
        'vx': lambda x: x['velocity'][:, 0],
        'vy': lambda x: x['velocity'][:, 1],
        'vz': lambda x: x['velocity'][:, -1],
    }
    if s in static:
        return static[s]

    if re.match('p[0-9]+', s):
        return lambda x: x['position'][:, int(s[1:])]
    if re.match('v[0-9]+', s):
        return lambda x: x['velocity'][:, int(s[1:])]

    try:
        a, b, *rest = re.split('(\d+)', s)
        if rest != ['']:
            warnings.warn(f'Ignored part of description: {rest}')
        return lambda x: x[a][:, int(b)]
    except ValueError:
        # Can't make sense of s, can only treat it as a direct accessor
        return lambda x: x[s]


def parse_dimension_description(s: str):
    if ',' in s:
        l, r = s.split(',')
        l = parse_single_dimension_description(l)
        r = parse_single_dimension_description(r)
        return lambda x: (l(x), r(x))
    else:
        return parse_single_dimension_description(s)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
