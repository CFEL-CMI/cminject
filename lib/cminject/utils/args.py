#!/usr/bin/env python3
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMInject
#
# Copyright (C) 2018,2020 CFEL Controlled Molecule Imaging group
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# If you use this program for scientific work, you should correctly reference it; see the LICENSE.md file for details.
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
import logging

import pyparsing
import re
import warnings
from typing import Union, Callable, Any, Iterable

import numpy as np

from cminject.utils.distributions import Distribution, parse_distribution

DimensionDescription = Union[str, Callable[[np.array], Any]]


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


def distribution_description(x):
    """
    Defines a custom argparse type for a distribution description.
    :param x: The string value to parse.
    :return: A :class:`Distribution` instance.
    """
    try:
        return parse_distribution(x)
    except pyparsing.ParseException as e:
        logging.info("PARSER ERROR: " + str(e))
        raise ValueError(f"{x} could not be parsed as a distribution! Enable --loglevel info to see the parser error.")


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


def parse_dimension_description(s: str) -> DimensionDescription:
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
        a, b, *rest = re.split(r'(\d+)', s)
        if rest != ['']:
            warnings.warn(f'Ignored part of description: {rest}')
        return lambda x: x[a][:, int(b)]
    except ValueError:
        # Can't make sense of s, can only treat it as a direct accessor
        return lambda x: x[s]


def parse_multiple_dimensions_description(s: str, n: Union[int, Iterable[int]] = None) -> (DimensionDescription, int):
    """
    Parses a dimension description from a string, which is either one name of a dimension like
      :func:`parse_dimension_description` accepts, or multiple of such names separated by commata.

    :param s: The string to parse as a dimension description.
    :param n: The number of expected/allowed components in the description.
      Either an integer or an iterable of integers. If None, no checks are performed.
    :return: A Callable that will return the appropriate component(s) from a given np.array. If multiple components were
      given as 's', this callable will return a corresponding k-tuple of components.
    """
    if ',' in s:
        components = s.split(',')
        if n is not None:
            assert (len(components) == n if type(n) is int else len(components) in n),\
                f"Dimension description of {n} components expected, reading {len(components)} components in '{s}'."

        fns = list(map(parse_dimension_description, components))
        # return a function that applies each function to the input, and returns a list of all of those applications
        return (lambda x: list(map(lambda fn: fn(x), fns))), len(components)
    else:
        assert (n == 1 if type(n) is int else 1 in n),\
            f"{n} components expected, but the dimension description '{s}' contains only one."
        return parse_dimension_description(s), 1
