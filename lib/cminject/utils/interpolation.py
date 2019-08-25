from typing import List, Tuple

import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d

from cminject.utils.cython_interpolation import interp2D, interp3D


class Interp2D:
    """
    A fast linear interpolator on a regular 2D grid with interpolation implemented in Cython,
    with a construction and call API like scipy.interpolate.RegularGridInterpolator. Wraps cython_interpolation.pyx.
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

    def __call__(self, t):
        # Get the size of the data grid in the output, X, and Y dimension
        D, X, Y = self.data.shape
        # Map the x,y position to the scale of the array indices
        x = (t[0] - self.min_x) / self.delta_x
        y = (t[1] - self.min_y) / self.delta_y

        return interp2D(self.data, x, y, D, X, Y)


class Interp3D:
    """
    Like Interp2D, but on a 3D grid.
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

    def __call__(self, t):
        # Get the size of the data grid in the output, X, Y, and Z dimension
        D, X, Y, Z = self.data.shape
        # Map the x,y,z position to the scale of the array indices
        x = (t[0] - self.min_x) / self.delta_x
        y = (t[1] - self.min_y) / self.delta_y
        z = (t[2] - self.min_z) / self.delta_z

        return interp3D(self.data, x, y, z, D, X, Y, Z)


def get_regular_grid_interpolator(grid, data):
    """
    Returns an appropriate regular grid interpolator instance based on the dimensionality. This is either an Interp2D,
        a Interp3D, or a RegularGridInterpolator instance, and Interp2D/Interp3D are much faster for 2D or 3D.

    The returned objects follow mostly the same API as scipy's RegularGridInterpolator, with the restriction that
        the coordinate argument must be a single coordinate and not a vector of coordinates.

    :param grid: The points defining the regular grid in n dimensions. Tuple/list of ndarray of float, with shapes
        (m1,), ..., (mn,)
    :param data: The data on the regular grid in n dimensions. Shape (m1, ..., mn, D). Interpolated data vectors with
        shape (D,) will be returned by the interpolator.
    :return: An interpolator instance that can be called like a function.
    """
    if len(grid) == 3:
        return Interp3D(grid, data)
    elif len(grid) == 2:
        return Interp2D(grid, data)
    else:
        return RegularGridInterpolator(grid, data)


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


def reconstruct_detector_hits(trajectories: np.array, z: float, interpolation_kind: str = 'linear') -> np.array:
    """
    Reconstructs measured quantities (e.g. x/y positions) at a given z position from list of trajectories,
    based on interpolating on each function piece of all trajectory curves.

    :param trajectories: The set of trajectories, which should be a (d, n)-shaped np.array, where d is the number
        of interpolated quantities and n is the number of recorded points in the trajectory.
    :param z: The z position to reconstruct the quantities at.
    :param interpolation_kind: The kind of interpolation to use. 'linear' by default.
        Refer to scipy.interpolate.interp1d for more information, as this is the interpolator used.
    :return: A (d, N)-shaped array of all interpolated quantities at the given Z position, where d is the number
        of interpolated quantities and 0 <= N is the number of points that could be reconstructed. N will be smaller
        than n if at least one trajectory did not pass through the given z position, and it can be larger than n if
        some trajectories passed through the given z position multiple times.
    """
    hit_xys = []

    for i, traj in enumerate(trajectories):
        zs = traj[-1]
        xys = traj[:-1]

        if xys.shape[1] < 2:  # no trajectory recorded
            continue

        z_split, split_idxs = split_at_inflections(zs)
        xys_split = np.split(xys, split_idxs, axis=1)

        for i, z_piece in enumerate(z_split):
            if z_piece.shape[0] < 2:  # skip over pieces that only contain 1 point
                continue

            z_piece, indices = np.unique(z_piece, return_index=True)
            xy_piece = xys_split[i][:, indices]

            ip = interp1d(z_piece, xy_piece, kind=interpolation_kind)
            try:
                hit_xys.append(ip(z))
            except ValueError:
                pass  # z was not in the function piece we're looking at

    return np.array(hit_xys).transpose()


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
