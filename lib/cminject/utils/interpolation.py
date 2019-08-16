import numpy as np
from scipy.interpolate import RegularGridInterpolator

from cminject.utils.cython_interpolation import interp2D, interp3D


class Interp2D:
    """
    A fast linear interpolator on a regular 2D grid with interpolation implemented in Cython,
    with a construction and call API like scipy.interpolate.RegularGridInterpolator. Wraps cython_interpolation.pyx.
    """
    def __init__(self, grid, data):
        x, y = grid
        self.data = np.moveaxis(data, -1, 0).copy(order='C')
        self.min_x, self.max_x = x[0], x[-1]
        self.min_y, self.max_y = y[0], y[-1]
        self.delta_x = (self.max_x - self.min_x) / (x.shape[0] - 1)
        self.delta_y = (self.max_y - self.min_y) / (y.shape[0] - 1)

    def __call__(self, t):
        D, X, Y = self.data.shape
        x = (t[0] - self.min_x) / self.delta_x
        y = (t[1] - self.min_y) / self.delta_y

        return interp2D(self.data, x, y, D, X, Y)


class Interp3D:
    """
    Like Interp2D, but on a 3D grid.
    """
    def __init__(self, c, data):
        x, y, z = c
        self.data = np.moveaxis(data, -1, 0).copy(order='C')
        self.min_x, self.max_x = x[0], x[-1]
        self.min_y, self.max_y = y[0], y[-1]
        self.min_z, self.max_z = z[0], z[-1]
        self.delta_x = (self.max_x - self.min_x) / (x.shape[0] - 1)
        self.delta_y = (self.max_y - self.min_y) / (y.shape[0] - 1)
        self.delta_z = (self.max_z - self.min_z) / (z.shape[0] - 1)

    def __call__(self, t):
        D, X, Y, Z = self.data.shape
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


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
