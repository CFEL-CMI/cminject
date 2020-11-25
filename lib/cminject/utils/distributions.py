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
- A collection of classes that serve as samplers of random distributions.
- Utility functions for parsing human-readable strings representing distributions, into instances of those classes.
"""

import abc
import pyparsing as pp
from typing import List, Union

import numpy as np


class Distribution(abc.ABC):
    """
    An abstract random 1D distribution that can generate and return a number of samples.
    Subclasses must implement the :meth:`generate` method to be instantiable.
    """
    @abc.abstractmethod
    def generate(self, n: int) -> np.array:
        """Return n samples of this distribution as a NumPy array."""
        pass

    def __add__(self, other: 'Distribution'):
        assert isinstance(other, Distribution), "Right operand must be Distribution instance"
        return _ArithDist(left=self, right=other, op='+')

    def __mul__(self, other: 'Distribution'):
        assert isinstance(other, Distribution), "Right operand must be Distribution instance"
        return _ArithDist(left=self, right=other, op='*')

    def __sub__(self, other: 'Distribution'):
        assert isinstance(other, Distribution), "Right operand must be Distribution instance"
        return _ArithDist(left=self, right=other, op='-')

    def __truediv__(self, other: 'Distribution'):
        assert isinstance(other, Distribution), "Right operand must be Distribution instance"
        return _ArithDist(left=self, right=other, op='/')

    def __repr__(self):
        return '<' + self.__str__() + '>'


class _ArithDist(Distribution):
    def __init__(self, left: 'Distribution', right: 'Distribution', op: str):
        assert op in ['+', '*', '-', '/'], f"Unknown operation '{op}'"
        self.left = left
        self.right = right
        self.op = op

    def generate(self, n: int) -> np.array:
        if self.op == '+':
            return self.left.generate(n) + self.right.generate(n)
        elif self.op == '*':
            return self.left.generate(n) * self.right.generate(n)
        elif self.op == '-':
            return self.left.generate(n) - self.right.generate(n)
        elif self.op == '/':
            return self.left.generate(n) / self.right.generate(n)

    def __str__(self):
        return f'({self.left.__str__()} {self.op} {self.right.__str__()})'


class GaussianDistribution(Distribution):
    """A Gaussian distribution with a fixed mean (mu) and standard deviation (sigma). See np.random.normal."""
    def __init__(self, mu: float, sigma: float):
        self.mu, self.sigma = mu, sigma

    def generate(self, n: int) -> np.array:
        return np.random.normal(self.mu, self.sigma, n)

    def __str__(self):
        return f'Gaussian(mu={self.mu},sigma={self.sigma})'


class UniformDistribution(Distribution):
    """A uniform distribution between a minimum and a maximum value. See np.random.uniform."""
    def __init__(self, min: float, max: float):
        self.min, self.max = min, max

    def generate(self, n: int) -> np.array:
        return np.random.uniform(self.min, self.max, n)

    def __str__(self):
        return f'Uniform(min={self.min},max={self.max})'


class LinearDistribution(Distribution):
    """A linear "distribution" between a minimum and a maximum value, i.e. generates a np.linspace."""
    def __init__(self, min: float, max: float):
        self.min, self.max = min, max

    def generate(self, n: int) -> np.array:
        return np.linspace(self.min, self.max, n)

    def __str__(self):
        return f'Linear(min={self.min},max={self.max})'


class DiracDeltaDistribution(Distribution):
    """
    The dirac delta distribution delta(x - value), samples of which are always the same constant value.
    """
    def __init__(self, value):
        """
        :param value: The constant value to return.
        """
        self.value = value

    def generate(self, n: int) -> np.array:
        return np.full((n,), self.value)

    def __str__(self):
        return f'DiracDelta({self.value})'


class NormOfDistributions(Distribution):
    """
    Returns the distribution of the p-norm of samples generated by n multiple other distributions.
    Can be used to, for example, generate the r distribution given the x and y distributions.

    Akin to a generalised chi distribution that is not necessarily based on normally distributed variables.
    """
    def __init__(self, *distributions: Distribution, p: int = 2):
        """
        :param distributions: The list of distributions to generate the norm distribution of.
          n (as used in the class docstring) is then equal to len(distributions).
        :param p: (Optional) Which p-norm to use. Defaults to 2 (i.e. the Euclidean norm).
        """
        self.distributions = distributions
        self.p = p

    def generate(self, n):
        dist_values = [np.abs(distribution.generate(n))**self.p for distribution in self.distributions]
        return np.sum(dist_values, axis=0)**(1/self.p)

    def __str__(self):
        return f'Norm_{self.p}({", ".join(str(d) for d in self.distributions)})'


# String -> Distribution parsing

def _tokens_to_dist(tokens):
    kind = tokens[0]
    assert kind in 'GLU', f'Only G/L/U are allowed as distribution kinds, you passed {kind}'
    if kind == 'G':
        assert len(tokens) == 6
        return GaussianDistribution(tokens[2], tokens[4])
    elif kind == 'L':
        assert len(tokens) == 6
        return LinearDistribution(tokens[2], tokens[4])
    elif kind == 'U':
        assert len(tokens) == 6
        return UniformDistribution(tokens[2], tokens[4])


def _tokens_to_dirac(tokens):
    return DiracDeltaDistribution(_tokens_to_float(tokens))


def _tokens_to_dist_norm(tokens):
    return NormOfDistributions(*tokens[0])


def _tokens_to_float(tokens):
    return float(tokens[0])


_fnumber_regex = r"[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?"
_fnumber_arg = pp.Regex(_fnumber_regex)\
    .setName('float')\
    .setParseAction(_tokens_to_float)
_const_dist = pp.Regex(_fnumber_regex)\
    .setName('distribution-constant')\
    .setParseAction(_tokens_to_dirac)

_plain_dist = (pp.Word('GLU') + '[' + _fnumber_arg + ',' + _fnumber_arg + ']')\
    .setName('distribution-simple')\
    .setParseAction(_tokens_to_dist)
_norm_dist = pp.nestedExpr('norm[', ']', content=pp.delimitedList(_plain_dist, ','))\
    .setName('norm-of-distributions')\
    .setParseAction(_tokens_to_dist_norm)
_dist = (_norm_dist | _plain_dist | _const_dist).setName('distribution')


def parse_distribution(s: str) -> Distribution:
    """
    Parse a string that describes a distribution and return a Distribution object. The allowed formats are given on the
    left, the corresponding Distribution subclass on the right.

    - G(mu,sigma) -> :class:`GaussianDistribution`
    - U(min,max)  -> :class:`UniformDistribution`
    - L(min,max)  -> :class:`LinearDistribution`

    Additionally, the following format is allowed:

    - norm(X1[,X2[,X3[,...]]])) -> :class:`NormOfDistributions`, where all Xi must match any of the formats specified
      here.

    :param s: The string to parse.
    :return: The constructed Distribution object, if parsing was successful. Otherwise raises an error.
    :raises: :class:`pp.ParseException`
    """
    result = _dist.parseString(s, parseAll=True)[0]
    return result


constant = DiracDeltaDistribution  # An intuitive and short alias
