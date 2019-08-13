import numpy as np
from scipy.interpolate import RegularGridInterpolator

from cminject.utils.cython_interpolation import interp2D, interp3D


class Interp2D:
    def __init__(self, c, v):
        x, y = c
        self.v = np.transpose(v, axes=(2, 0, 1)).copy(order='C')
        self.min_x, self.max_x = x[0], x[-1]
        self.min_y, self.max_y = y[0], y[-1]
        self.delta_x = (self.max_x - self.min_x) / (x.shape[0] - 1)
        self.delta_y = (self.max_y - self.min_y) / (y.shape[0] - 1)

    def __call__(self, t):
        V, X, Y = self.v.shape
        x = (t[0] - self.min_x) / self.delta_x
        y = (t[1] - self.min_y) / self.delta_y

        return interp2D(self.v, x, y, V, X, Y)


class Interp3D:
    def __init__(self, c, v):
        x, y, z = c
        self.v = np.transpose(v, axes=(3, 0, 1, 2)).copy(order='C')
        self.min_x, self.max_x = x[0], x[-1]
        self.min_y, self.max_y = y[0], y[-1]
        self.min_z, self.max_z = z[0], z[-1]
        self.delta_x = (self.max_x - self.min_x) / (x.shape[0] - 1)
        self.delta_y = (self.max_y - self.min_y) / (y.shape[0] - 1)
        self.delta_z = (self.max_z - self.min_z) / (z.shape[0] - 1)

    def __call__(self, t):
        V, X, Y, Z = self.v.shape
        x = (t[0] - self.min_x) / self.delta_x
        y = (t[1] - self.min_y) / self.delta_y
        z = (t[2] - self.min_z) / self.delta_z

        return interp3D(self.v, x, y, z, V, X, Y, Z)


def get_interpolator(c, v):
    if len(c) == 3:
        return Interp3D(c, v)
    elif len(c) == 2:
        return Interp2D(c, v)
    else:
        return RegularGridInterpolator(c, v)