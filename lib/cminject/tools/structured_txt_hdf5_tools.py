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
import logging
import re
import warnings
from collections import OrderedDict
from typing import Tuple, List, Dict, Union

import h5sparse
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

gte2spaces = re.compile(r'\s{2,}')

def _mirror_around_axis(arr, axis=0, flipsign=False):
    """
    Mirrors a 2D array around an axis and returns the resulting array.
    The "middle" column in the resulting array is not duplicated,
    i.e. the axis is considered to cut through that column.

    :param arr: The array. Must be 2-dimensional, i.e. len(arr.shape) must be 2.
    :param axis: The axis to flip along. Must be 0 or 1.
    :param flipsign: Whether to flip the sign of the values in the array on the mirrored side.
    :return: An (2m-1, n)-shaped (for axis=0) or a (m, 2n-1)-shaped array (for axis=1),
        where the mirrored part is constructed as described in the docstring above.
    """
    dimensions = len(arr.shape)
    if dimensions not in [1, 2]:
        raise ValueError(f"Cannot mirror a {arr.shape}-shaped array with this method.")
    elif dimensions == 2:
        if axis not in [0, 1]:
            raise ValueError("Can only mirror around axis 0 or 1 in a 2-D array!")
    elif dimensions == 1 and axis != 0:
        raise ValueError("Can only mirror around axis 0 in a 1-D array!")

    # Mirror around the axis and optionally flip the sign of the values
    mir = np.flip(arr, axis=axis)
    if flipsign:
        mir = -mir  # velocity in the direction of the mirrored dimension must be mirrored

    if dimensions == 2:  # Case for 2D
        m, n = arr.shape
        stacked = np.stack([mir, arr], axis=axis).reshape((m * 2, n))  # Stack the arrays along the mirror axis
        sliced = stacked[list(range(m)) + list(range(m + 1, 2 * m)), :]  # Cut out the twice-occurring middle column
        return sliced
    else:  # Case for 1D
        m = arr.shape[0]
        stacked = np.concatenate([mir, arr])  # Stack the arrays along the mirror axis
        sliced = stacked[list(range(m)) + list(range(m + 1, 2 * m))]  # Cut out the twice-occurring middle value
        return sliced


def txt_to_hdf5(infile_name: str, outfile_name: str, dimensions: int = 3, mirror: bool = False,
                attributes: Union[Dict, None] = None) -> None:
    """
    A function that reads a structured .txt file and stores a sparse specially constructed HDF5 file,
    while keeping metadata from the original file. The general structure of the .txt input file must be as follows::

        % <Some meta information>
        % <Some more meta information>
        % X  Y  Z  A (m)  B (s)  C ...
        0.0  0.0  0.0  1.0  2.0  3.0 ...
        0.1  0.0  0.0  4.0  5.0  6.0 ...
        ...

    This matches the .txt output format of COMSOL for 3D grids.

    All but the last line starting with % at the beginning of the file are considered to be a general textual
    description of the file (metadata about the file itself), and are stored.

    The last line starting with % is considered the header of the columns, where each column name has to be separated
    by at least two spaces from the other column names. A single space is considered to still be part of the same
    column, and the second part of the column name is - if present - considered to be the unit of the column.
    This will be stored as metadata too.

    The number of columns after X/Y/Z is arbitrary.

    The constructed HDF5 file will have the following entries:

    - index: Each of these is the set of indices found for each dimension, in order of occurrence

      - x
      - y
      - z, if present in the original file

    - data: Each of these is one full row of values, in X*Y*Z order (X and Z are flipped wrt. the original txt file!)

      - <one for each column, matching the column name from the file>

    The metadata (the lines starting with %) is stored as an attribute on the HDF5 root.
    Columns specifying a unit (separated by one space) have this information attached as metadata on data/<...>.

    :param infile_name: The input file name. Must be a .txt file matching the format described above, otherwise
        the behavior will be undefined.
    :param outfile_name: The output file name.
    :param dimensions: The number of spatial dimensions the txt file is defined in. 3 by default, will be 2 or 3 for
        most cases.
    :param mirror: A flag for symmetry around an axis: If true, the entire field will be duplicated and mirrored around
        an axis in the first dimension that's positioned the minimal position in the first dimension.
    :param attributes: Optionally, a dictionary of attributes to attach to the output HDF5 file (as HDF5 attributes on
        the root node). Useful to store global properties of the field or anything else that tools or people who later
        read this file might need.
    :return: None.

    """
    f = open(infile_name, 'r')
    # strictly speaking we only need a set for storing indices, but there is no OrderedSet in base python3
    index = [OrderedDict() for _ in range(dimensions)]
    values = {}

    # Handle files with either headers or no headers
    has_headers = f.read(1) == '%'  # The file is assumed to only have headers if it begins with a % sign
    f.seek(0)

    if has_headers:
        metadata, headers = [], []
        prev_line = ""
        last_pos = 0

        line = f.readline()
        while line:
            if line.startswith('%'):
                metadata.append(line)
                prev_line = line
            elif not line.startswith('%') and prev_line.startswith('%'):
                # The previous line is the "header line" if this one does not start with a %
                headers = [tuple(x.split()) for x in re.split(gte2spaces, prev_line[1:].strip())]
                # Discard the headers for the spatial dimensions
                headers = headers[dimensions:]
                # Store an empty array in values for each header name
                for header in headers:
                    values[header[0]] = []
                break
            last_pos = f.tell()
            line = f.readline()

        f.seek(last_pos)
    else:
        warnings.warn(
            "WARNING: The .txt file to convert has no headers, and so we can't assign names for columns. "
            "If you know the relevant info, it's best to add these headers by hand and run the conversion again.")
        metadata = []
        # Read the first line and rewind again
        first_line = f.readline()
        f.seek(0)
        # Get the number of columns by splitting the line (on whitespace), and construct the header as just 0-based ints
        headers = [(i,) for i in range(len(first_line.split()) - dimensions)]
        for header in headers:
            values[header[0]] = []

    # Needed below, to recognise whether we need to flip the index or not
    first_cols = None
    flip_indices = None
    guessed_flip_indices_flag = False

    # From here on we assume every line to only contain data
    for line_no, line in enumerate(f):
        # After the header handling: Split current line on spaces, update index and values
        cols = line.split()
        for d in range(dimensions):
            # Add an entry to the index for the current row, for each dimension
            index[d][float(cols[d])] = None  # We only use the key, value can be None
        for i in range(len(headers)):
            values[headers[i][0]].append(float(cols[dimensions+i]))

        # To recognise whether we need to flip the index order or not, based on the first two data lines:
        if first_cols is None:
            # Remember the first set of columns, this will only occur once
            first_cols = cols
        elif not guessed_flip_indices_flag:
            # This will also only occur once, and set the "flip_indices" flag based on whether the first column has
            # changed between the first and second lines or not. If it has, we need to flip the index. Otherwise not.
            flip_indices = float(first_cols[0]) != float(cols[0])
            if flip_indices:
                logging.info(
                    "NOTE: Recognised that conversion needs to flip the index. You can read up on what this means in "
                    "the code, but most likely this will work exactly as intended for you.")
            guessed_flip_indices_flag = True

    # Construct numpy arrays for each collected index set per dimension
    index = [np.array(list(i.keys()), dtype=float) for i in index]
    index_n = [len(i) for i in index]

    for i in range(len(headers)):
        val_arr = np.array(values[headers[i][0]], dtype=np.float)
        if flip_indices:
            # Transpose the whole array to match "normal" X/Y/Z/... iteration order. The order that the dimensions
            # increment in is inverted in the txt file format, so transpose it here.
            val_arr = val_arr.reshape(tuple(reversed(index_n))).transpose().reshape((np.prod(index_n),))

        if mirror:
            # Determine if we need to flip the sign on the mirrored side (only for v_r)
            flipsign = headers[i][0] in ['v_r', 'V_R', 'u', 0]  # TODO better conditions to check for?
            # Mirror the value array along that axis, possibly flipping the mirrored values' signs
            # TODO what if flip_indices *and* mirror is True? is this correct then? I think so but I'm unsure  - Simon
            val_arr = _mirror_around_axis(val_arr.reshape(index_n), axis=0, flipsign=flipsign).flatten()

        # Construct sparse matrix and store it in values
        val_mat = csr_matrix(np.nan_to_num(val_arr))
        values[headers[i][0]] = val_mat

    if mirror:
        # Update the index array for the first dimension
        index[0] = _mirror_around_axis(index[0], axis=0, flipsign=True)

    _save_to_hdf5(attributes, dimensions, headers, index, metadata, outfile_name, values)

def txt_to_hdf5_stark(field_name: str, grad_name: str, outfile_name: str, dimensions: int = 2, mirror: bool = False):
    """
    A funciton that reads field and gradient txt. files describing stark deflectors. This functions is needed because Comsol text files
    for stark deffectors are different from txt files expected in txt_to_hdf5. The structure of the txt file expected:
        Initial y-coordinate [mm]    initial x-coordinate [mm]    initial z-coordinate [mm]
        stepsize_y           [mm]    stepsize_x           [mm]    stepsize_z           [mm]
        #grid_points in y-direction  #grid_points in x-direction  #grid_points in z-direction
        f(x0,y0) [v/m]
        f(x1,y0) [v/m]
        .
        .
        f(x0,y1)
        f(x1, y1)
        .
        .
        f(xn, yn)

    :param field_name: the filename of the field strength
    :param grad_name: the filename of the gradient
    :param dimensions: The number of spatial dimensions of the grid. 2 by default, will be 2 or 3 for
        most cases.
    :param mirror: A flag for symmetry around an axis: If true, the entire field will be duplicated and mirrored around
        an axis in the first dimension that's positioned the minimal position in the first dimension. Not yet applied in this script.
    :return: None.

    The constructed HDF5 file will have the following entries:

    - index: Each of these is the set of indices found for each dimension, in order of occurrence

      - x
      - y
      - z, if present in the original file

    - data: Each of these is one full row of values, in X*Y*Z order (X and Z are flipped wrt. the original txt file!)
    """
    # reading header (it should be the same for the field and the gradient files)
    info = np.genfromtxt(field_name, dtype = float, max_rows = 3)
    field = np.genfromtxt(field_name, dtype = float, skip_header = 3)
    gradient = np.genfromtxt(grad_name, dtype = float, skip_header = 3)

    # creating the grids
    bins_y = int(info[2,0]) # number of points in the y-direction
    bins_x = int(info[2,1]) # number of points in the x-direction
    y0 = info[0,0] # the starting point in the y-direction
    x0 = info[0,1] # the starting point in the x-direction
    s_y = info[1,0] # stepsize in the y-direction
    s_x = info[1,1] # stepsize in the x-direction
    index_x = np.linspace(x0, x0 + s_x*bins_x, num = bins_x)
    index_y = np.linspace(y0, y0 + s_y*bins_y, num = bins_y)
    # if the field is only one dimensional it will be stored in an array of shape (n,). Hence, we must change it
    try:
        np.shape(field)[1]
    except:
        field = field.reshape(-1,1)

    # Create names for the coloumns of the data we have
    headers_g = ['ex', 'ey', 'ez']
    column_index_g = [0, 1, 2]
    if np.shape(field)[1] == 3:
        headers_f = ['fx', 'fy', 'fz']
        column_index_f = [3,4,5]
    else:
        headers_f = ['f_norm']
        column_index_f = [3]
    """ Flipping the matrix of the field and gradient to be
    in an increasing order with respect to the y variable
    Note: I need to edit this in case we have a z-dimension, i.e. fringe field
    """
    # creating an indicator matrix
    a = np.arange(0, bins_x*bins_y)
    a = np.reshape(a, (bins_x, bins_y))
    indicator = np.zeros((bins_x*bins_y,), dtype = int)
    cnt = 0
    for x in np.nditer(a, order = 'F'):
        indicator[cnt] = x
        cnt += 1
    # saving the data to an hdf5 file
    with h5sparse.File(outfile_name, 'w') as hp:
        hp.create_dataset('index/0', data = index_x)
        hp.create_dataset('index/1', data = index_y)
        group = hp.create_group('data')
        for f in range(0, np.shape(field)[1]):
            #group1 = hp.create_group('data/' + )
            data = field[:,f]
            flipped = np.zeros(np.shape(data), dtype = float)
            flipped[indicator[0:]] = data[0:]
            dat1 = group.create_dataset(headers_f[f], data = flipped)
            dat1.attrs['unit'] = 'v/m'
            dat1.attrs['column_index'] = column_index_f[f]

        for g in range(0, np.shape(gradient)[1]):
            #group2 = hp.create_group('data/' + )
            data = gradient[:, g]

            flipped = np.zeros(np.shape(data), dtype = float)
            flipped[indicator[0:]] = data[0:]
            dat2 = group.create_dataset(headers_g[g], data = flipped)
            dat2.attrs['unit'] = 'v/m^2'
            dat2.attrs['column_index'] = column_index_g[g]

        """
        code to flip the way the field is saved
        val_arr = np.array(values[headers[i][0]], dtype=np.float)
        if flip_indices:
            # Transpose the whole array to match "normal" X/Y/Z/... iteration order. The order that the dimensions
            # increment in is inverted in the txt file format, so transpose it here.
            val_arr = val_arr.reshape(tuple(reversed(index_n))).transpose().reshape((np.prod(index_n),))
        """

def _save_to_hdf5(attributes, dimensions, headers, index, metadata, outfile_name, values):
    # Use h5sparse to save (lots of) space when storing (lots of) NaN or 0.0 values
    with h5sparse.File(outfile_name) as h5f:
        # Store the metadata from the .txt file and the passed attributes as attributes on the root of the hdf5 file
        h5f.attrs['metadata'] = "\n".join(metadata)
        if attributes:
            for k, v in attributes.items():
                h5f.attrs[k] = v

        # Create a dataset for each column, store the np.array and the original index of the column
        for i in range(len(headers)):
            dataset = h5f.create_dataset(f'data/{headers[i][0]}', data=values[headers[i][0]])
            dataset.attrs['column_index'] = i  # store this to preserve column order
            if len(headers[i]) > 1:
                # If the header for this column has more than one part, it's assumed the second part is the unit
                # and this is also stored as metadata
                dataset.attrs['unit'] = headers[i][1]

        # Store the index definitions for each spatial dimension
        for i in range(dimensions):
            h5f[f'index/{i}'] = index[i]


def hdf5_to_data_frame(filename: str) -> pd.DataFrame:
    """
    Reads in an HDF5 file (as written by txt_to_hdf5) and returns a pandas DataFrame matching it.
    Note that not all metadata stored in the HDF5 file is also attached to the DataFrame as of now.

    :param filename: The HDF5 file to read in. Must be written by or adhere to the format that txt_to_hdf5 defines.
    :return: A pandas DataFrame matching the data stored in the HDF5 file.
    """
    with h5sparse.File(filename, 'r') as h5f:
        dimensions = len(h5f['index'])
        index = [h5f[f'index/{i}'] for i in range(dimensions)]

        column_indices = []
        column_headers = []
        column_values = []

        # Populate the values dict from the columns in 'data'
        for column_name in h5f['data']:
            col = h5f['data'][column_name]
            column_headers.append(column_name)
            column_indices.append(col.attrs['column_index'])

            values = col[()]
            if not isinstance(values, np.ndarray):
                values = values.toarray()
            column_values.append(values)

        # Reorder the column data and headers according to the (inverse) permutation gathered from the metadata
        column_indices = np.argsort(column_indices)  # argsort constructs the inverse permutation
        column_headers = np.array(column_headers, dtype=str)[column_indices]
        # Values needs to be (n_d x n_c), not (n_c x n_d) for the DataFrame constructor
        # where n_d is the number of data rows and n_c is the number of columns. So swap the axes.
        column_values = np.vstack(column_values)[column_indices].swapaxes(0, 1)

        # Generate a MultiIndex from the product of the indices for each spatial dimension
        # TODO index names should be reconstructed from the file
        idx = pd.MultiIndex.from_product(index, names=['x', 'y', 'z'] if dimensions == 3 else ['r', 'z'])
        # Construct the DataFrame from the values with the index
        df = pd.DataFrame(data=column_values, columns=column_headers, index=idx)

        return df


def data_frame_to_data_grid(df: pd.DataFrame) -> Tuple[List[np.array], np.array]:
    """
    Turns a DataFrame (as constructed by hdf5_to_data_frame) into a data grid that can be used with the
    scipy.interpolate.RegularGridInterpolator.
    The result can be passed on like: RegularGridInterpolator((x,y,z), data_grid).

    :param df: The pandas DataFrame. Must be constructed by or adhere to the format that hdf5_to_data_frame defines.
    :return: A 2-tuple like ([x, y, ...], data_grid), where the first component is a list of all index arrays,
        and the second component is the data grid matching this index.
    """
    # Gather index and number of indices along each axis from the DataFrame's MultiIndex
    index = [level.values for level in df.index.levels]
    index_n = [len(df.index.levels[i]) for i in range(len(df.index.levels))]

    # Create the data grid and shape it into expected n-dimensional form
    new_shape = tuple(index_n + [-1])
    data_grid = np.nan_to_num(df.to_numpy().reshape(new_shape))
    # Return the index array and the data grid. Both are needed to construct a RegularGridInterpolator.
    return index, data_grid


def hdf5_to_data_grid(filename: str) -> Tuple[List[np.array], np.array]:
    """
    Shortcut function to read in an HDF5 file (as written by txt_to_hdf5) and turn it into a data grid.
    Refer to hdf5_to_data_frame and data_frame_to_data_grid for further info; this function is defined purely
    in terms of the composition of the two.

    :param filename: The HDF5 file to read in.
    :return: A 4-tuple of numpy arrays: (x, y, z, data_grid).
    """
    return data_frame_to_data_grid(hdf5_to_data_frame(filename))

if __name__ == '__main__':

    e_file = '../../../../gt.grad.txt'
    f_file = '../../../../Ec.norm.txt'
    txt_to_hdf5_stark(f_file, e_file, 'test')
    #d = {'J': 2, 'Ka': 2, 'Kc': 2, 'M': 2, 'Isomer':0}
    #pos = np.array([1,2,3])
    #p = Molecule(34., d)
    #f = B_StarkField(f_file, e_file)
    #f.energy_interpolate(p)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
