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

import logging
from typing import List, Tuple

import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d

from cminject.utils.cython_interpolation import interp2D, interp3D


# --- For n-linear interpolation on a regular grid

class Interp2D:
    """
    A fast linear interpolator on a regular 2D grid with interpolation implemented in Cython,
    with a construction and call API like scipy.interpolate.RegularGridInterpolator, except
    that it accepts a (2,)-shaped np.array directly at call and not a 2-tuple of (x,y) coordinate
    arrays. Wraps cython_interpolation.pyx.
    """
    def __init__(self, grid, data):
        x, y = grid
        # Move the output dimension axis into the first spot and copy it to be C-contiguous, enabling the interpolating
        # code to easily find the slice of the array to interpolate on for each output dimension.
        self.data = np.moveaxis(data, -1, 0).copy(order='C')
        # Calculate and store the extent and spacing for X/Y.
        self.min_x, self.max_x = x[0], x[-1]
        self.min_y, self.max_y = y[0], y[-1]
        self.delta_x = (self.max_x - self.min_x) / (x.shape[0] - 1)
        self.delta_y = (self.max_y - self.min_y) / (y.shape[0] - 1)
        # A vector of appropriate output size to copy, modify and return from each call
        self._nanvec = np.full(self.data.shape[0], np.nan, dtype=np.float64)

    def __call__(self, t):
        """
        Return the interpolated data at the 3D coordinate coord.

        :param coord: The coordinate to interpolate at, as a (2,)-shaped np.array.
        :return: An (nd,)-shaped np.array with the interpolated results.
        """
        # Get the size of the data grid in the output, X, and Y dimension
        D, X, Y = self.data.shape
        # Map the x,y position to the scale of the array indices
        x = (t[0] - self.min_x) / self.delta_x
        y = (t[1] - self.min_y) / self.delta_y

        out = self._nanvec.copy()
        interp2D(self.data, x, y, D, X, Y, out)
        return out


class Interp3D:
    """
    Like Interp2D, but on a 3D regular grid.
    """
    def __init__(self, c, data):
        x, y, z = c
        # Move the output dimension axis into the first spot and copy it to be C-contiguous, enabling the interpolating
        # code to easily find the slice of the array to interpolate on for each output dimension.
        self.data = np.moveaxis(data, -1, 0).copy(order='C')
        # Calculate and store the extent and spacing for X/Y/Z.
        self.min_x, self.max_x = x[0], x[-1]
        self.min_y, self.max_y = y[0], y[-1]
        self.min_z, self.max_z = z[0], z[-1]
        self.delta_x = (self.max_x - self.min_x) / (x.shape[0] - 1)
        self.delta_y = (self.max_y - self.min_y) / (y.shape[0] - 1)
        self.delta_z = (self.max_z - self.min_z) / (z.shape[0] - 1)
        # A vector of appropriate output size to copy, modify and return from each call
        self._nanvec = np.full(self.data.shape[0], np.nan, dtype=np.float64)

    def __call__(self, coord: np.array):
        """
        Return the interpolated data at the 3D coordinate coord.

        :param coord: The coordinate to interpolate at, as a (3,)-shaped np.array.
        :return: An (nd,)-shaped np.array with the interpolated results.
        """
        # Get the size of the data grid in the output, X, Y, and Z dimension
        D, X, Y, Z = self.data.shape
        # Map the x,y,z position to the scale of the array indices
        x = (coord[0] - self.min_x) / self.delta_x
        y = (coord[1] - self.min_y) / self.delta_y
        z = (coord[2] - self.min_z) / self.delta_z

        out = self._nanvec.copy()
        interp3D(self.data, x, y, z, D, X, Y, Z, out)
        return out


class InterpND:
    """
    A simple wrapper for n-dimensional n-linear interpolation on a regular grid. Wraps
    scipy.interpolate.RegularGridInterpolator directly except for an additional `tuple()` call
    around the passed coordinate to interpolate at. This is to make the API coherent with Interp2D
    and Interp3D, which directly accept NumPy arrays, as this is faster for single-coordinate
    interpolation which Interp2D/Interp3D are optimized for.
    """
    def __init__(self, *args, **kwargs):
        self.interpolator = RegularGridInterpolator(*args, **kwargs)

    def __call__(self, t):
        return self.interpolator(tuple(t))


def get_regular_grid_interpolator(grid, data):
    """
    Returns an appropriate regular grid interpolator instance based on the dimensionality.
    This is either an Interp2D, an Interp3D, or a RegularGridInterpolator instance, and
    Interp2D/Interp3D are much faster for 2D or 3D.

    The returned objects follow mostly the same API as scipy's RegularGridInterpolator, with the
    restriction that the coordinate argument must be a single coordinate and not a vector of
    coordinates.

    :param grid: The points defining the regular grid in n dimensions. Tuple/list of ndarray of
        float, with shapes (m1,), ..., (mn,)
    :param data: The data on the regular grid in n dimensions. Shape (m1, ..., mn, D). Interpolated
        data vectors with shape (D,) will be returned by the interpolator.
    :return: An interpolator instance that can be called like a function.
    """
    if len(grid) == 3:
        return Interp3D(grid, data)
    elif len(grid) == 2:
        return Interp2D(grid, data)
    else:
        return InterpND(grid, data)


# --- For interpolation along particle trajectories

def split_at_inflections(a: np.array) -> Tuple[List[np.array], np.array]:
    """
    Splits a np.array into multiple pieces based on where the pairwise differences of the elements change sign,
    i.e. the inflection points of a curve going through the points. The resulting pieces are either monotonically
    increasing or decreasing, but not necessarily strictly (i.e. there are no splits at the same value occurring twice).

    Useful for interpolating through a curve since multiple interpolating functions for one curve can be generated
    from the result.

    :param a: The array to split as described.
    :return: A 2-tuple of:
        - A list of arrays, each being monotonically increasing or decreasing (not necessarily strictly).
        - The indices where the splits occurred.
    """
    # Get the pairwise differences
    diffs = np.sign(np.diff(a))
    # Find the positions of where the sign of those differences, also compared pairwise, changed
    changes = ((np.roll(diffs, 1) - diffs) != 0)
    changes[0] = False  # don't compare first to last
    # Get the indices of where the changes occurred
    idxs = np.where(changes)[0]  # [0] since there is only one dimension
    # Split the input array at those indices, and return the split array along with the indices
    return np.split(a, idxs), idxs


def reconstruct_detector_hits(trajectories: List[np.array], xdims: List[int], zdim: int,
                              zs: List[float], interpolation_kind: str = 'linear') -> List[np.array]:
    """
    Reconstructs measured quantities (e.g. x/y positions) at a given z position from list of trajectories,
    based on interpolating on each function piece of all trajectory curves.

    :param trajectories: A list of trajectories. Each value contained should should be a (d, n)-shaped np.array,
        where d is the number of interpolated quantities and n is the number of recorded points in the trajectory.
    :param xdims: The indices (in the trajectory array) of the quantities to reconstruct.
    :param zdim: The index (in the trajectory array) of the dimension corresponding to z.
    :param zs: The z positions to reconstruct the quantities at.
    :param interpolation_kind: The kind of interpolation to use. 'linear' by default.
        Refer to scipy.interpolate.interp1d for more information, as this is the interpolator used.
    :return: A list of (d, N_i)-shaped arrays of all interpolated quantities at each given z position,
        where d is the number of interpolated quantities and 0 <= N_i is the number of points that could
        be reconstructed at the i-th z position. N_i will be smaller than n if at least one trajectory did not pass
        through the i-th z position, and it can be larger than n if some trajectories passed through the i-th
        z position multiple times.
    """
    rec_zs = np.array(zs)
    hit_xys = []

    for i, traj in enumerate(trajectories):
        if (i % 1000) == 0:
            logging.info(f"...{i}...")

        zs = traj[zdim]
        xys = traj[xdims]

        if xys.shape[1] < 2:  # no trajectory recorded
            continue

        z_split, split_idxs = split_at_inflections(zs)
        xys_split = np.split(xys, split_idxs, axis=1)

        for j, z_piece in enumerate(z_split):
            if z_piece.shape[0] < 2:  # skip over pieces that only contain 1 point
                continue

            z_piece, indices = np.unique(z_piece, return_index=True)
            xy_piece = xys_split[j][:, indices]

            ip = interp1d(z_piece, xy_piece, kind=interpolation_kind, bounds_error=False, fill_value=np.nan)
            try:
                hit_xys.append(ip(rec_zs))
            except ValueError:
                pass  # z was not in the function piece we're looking at

    ds = []
    hit_xys = np.array(hit_xys)
    for j in range(hit_xys.shape[-1]):
        d = hit_xys[..., j]
        d = d[~np.isnan(d).any(axis=1)].transpose()
        ds.append(d)
    return ds

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
