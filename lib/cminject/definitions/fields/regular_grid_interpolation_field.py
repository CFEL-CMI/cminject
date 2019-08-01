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

from abc import ABC
from typing import Tuple
from functools import partial

import numpy as np
from cminject.definitions.base import Field, Particle
from cminject.tools.structured_txt_hdf5_tools import hdf5_to_data_grid
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import solve_banded
from .interpolation import fast_cubic_spline as _spline

def f(coordinate, rGrid, zGrid, vr_c, vz_c, p_c):
    Vr=_spline.interpolate_2d(coordinate[0], coordinate[1], min(rGrid), max(rGrid), min(zGrid), max(zGrid), vr_c)
    Vz=_spline.interpolate_2d(coordinate[0], coordinate[1], min(rGrid), max(rGrid), min(zGrid), max(zGrid), vz_c)
    p=_spline.interpolate_2d(coordinate[0], coordinate[1], min(rGrid), max(rGrid), min(zGrid), max(zGrid), p_c)
#    print([Vr,Vz,p])
    return [Vr,Vz,p]


class RegularGridInterpolationField(Field, ABC):
    """
    A generic base class for fields that calculate a force based on interpolation within some n-dimensional
    regular grid. This class provides an `interpolate` method; the actual calculation of acceleration has to be
    defined in a subclass and should use the result of this method.
    """
    def __init__(self, filename: str, interpolator_type: str = 'regular_grid'):
        """
        The constructor for RegularGridInterpolationField.

        :param filename: The filename of an HDF5 file in the format as constructed by `tools/txt_to_hdf5`.
        :param interpolator_type: The interpolator to use. Options are 'regular_grid' and 'cython_cubic_spline',
            'regular_grid' is used by default.
        """
        # Construct the interpolator from the HDF5 file passed by file name
        self.filename = filename
        data_index, data_grid = hdf5_to_data_grid(self.filename)
        self.number_of_dimensions = len(data_index)
        
        if interpolator_type == 'regular_grid':
            self.interpolator = RegularGridInterpolator(tuple(data_index), data_grid)
        elif interpolator_type == 'cython_cubic_spline':
            if self.number_of_dimensions!=2:
                raise ValueError(f"Number of dimensions has to be 2 for this interpolator: {interpolator_type}")
            self.interpolator = self.construct_spline_interpolator(data_index, data_grid)
        else:
            raise ValueError(f"Unknown interpolator type: {interpolator_type}")

        # Assuming that Z is always the last dimension, construct the Z boundary from it
        minima = list(np.min(d) for d in data_index)
        maxima = list(np.max(d) for d in data_index)
        self._z_boundary = (minima[self.number_of_dimensions - 1], maxima[self.number_of_dimensions - 1])
        self.grid_boundary = np.array(list(zip(minima, maxima)))
        super().__init__()
    
    @staticmethod
    def construct_spline_interpolator(data_index, data_grid):
        def cal_coefs(a, b, y, c=None, alpha=0, beta=0):
            '''
            Return spline coefficients 
        
            Parameters
            ----------
            a : float
                lower bound of the grid.
            b : float
                upper bound of the grid.
            y : ndarray
                actual function value at grid points.
            c : (y.shape[0] + 2, ) ndarray, optional
                ndarry to be written
            alpha : float
                Second-order derivative at a. Default is 0.
            beta : float
                Second-order derivative at b. Default is 0.
        
            Returns
            -------
            out : ndarray
                Array of coefficients.
            '''
            n = y.shape[0] - 1
            h = (b - a)/n
        
            if c is None:
                c = np.zeros((n + 3, ))
                ifreturn = True
            else:
                assert(c.shape[0] == n + 3)
                ifreturn = False
        
            c[1] = 1/6 * (y[0] - (alpha * h**2)/6)
            c[n + 1] = 1/6 * (y[n] - (beta * h**2)/6)
        
            # ab matrix here is just compressed banded matrix
            ab = np.ones((3, n - 1))
            ab[0, 0] = 0
            ab[1, :] = 4
            ab[-1, -1] = 0
        
            B = y[1:-1].copy()
            B[0] -= c[1]
            B[-1] -=  c[n + 1]
        
            c[2:-2] = solve_banded((1, 1), ab, B)
        
            c[0] = alpha * h**2/6 + 2 * c[1] - c[2]
            c[-1] = beta * h**2/6 + 2 * c[-2] - c[-3]
        
            if ifreturn:
                return(c)
        rGrid,zGrid=data_index
        vr_c_tmp=np.zeros((len(rGrid)+2,len(zGrid)))
        cal_coefs(min(rGrid), max(rGrid), data_grid[...,0], vr_c_tmp)
        vr_c=np.zeros((len(rGrid)+2,len(zGrid)+2))
        cal_coefs(min(zGrid), max(zGrid), vr_c_tmp.T, vr_c.T)
        vz_c_tmp=np.zeros((len(rGrid)+2,len(zGrid)))
        cal_coefs(min(rGrid), max(rGrid), data_grid[...,1], vz_c_tmp)
        vz_c=np.zeros((len(rGrid)+2,len(zGrid)+2))
        cal_coefs(min(zGrid), max(zGrid), vz_c_tmp.T, vz_c.T)
        p_c_tmp=np.zeros((len(rGrid)+2,len(zGrid)))
        cal_coefs(min(rGrid), max(rGrid), data_grid[...,2], p_c_tmp)
        p_c=np.zeros((len(rGrid)+2,len(zGrid)+2))
        cal_coefs(min(zGrid), max(zGrid), p_c_tmp.T, p_c.T)
        return partial(f, zGrid=zGrid, rGrid=rGrid, vr_c=vr_c, vz_c=vz_c, p_c=p_c)
            

    @property
    def z_boundary(self) -> Tuple[float, float]:
        return self._z_boundary

    def is_particle_inside(self, particle: Particle, time: float) -> bool:
        return np.all(particle.spatial_position >= self.grid_boundary[:, 0]) and\
               np.all(particle.spatial_position <= self.grid_boundary[:, 1])

    def set_number_of_dimensions(self, number_of_dimensions: int):
        if number_of_dimensions != self.number_of_dimensions:
            raise ValueError(
                f"This field was constructed from a {self.number_of_dimensions}-dimensional grid and is "
                f"thus incompatible with a {number_of_dimensions}-dimensional simulation!"
            )


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""
