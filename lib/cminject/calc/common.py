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
Calculation functions that are potentially useful for several different physics models.
"""

import math
import numpy as np
import numba

inf = np.inf
ninf = np.NINF


@numba.jit("float64[::1](float64[::1])", nopython=True)
def erf_vec(a):
    """
    A 1D vectorized version of the error function, implemented in terms of math.erf.
    Mainly exists so we can call it for a 1D vector from Numba, which lacks support for scipy.special.erf.

    :param a: A 1D array of float64, to calculate erf component-wise for.
    :return: The component-wise application of erf to a.
    """
    return np.array([math.erf(x) for x in a])


def is_finite(x: float):
    """
    Like np.isfinite for a single floating point value (much faster). Returns True if x is not equal to +inf,
    -inf and if x is not NaN.

    :param x: The number to check.
    :return: True if x is finite, False otherwise.
    """
    return x != inf and x != ninf and not (x != x)


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
