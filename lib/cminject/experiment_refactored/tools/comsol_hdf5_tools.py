"""
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
"""
from typing import Tuple

import h5py
import h5sparse

from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd

import re

from collections import OrderedDict

gte2spaces = re.compile(r'\s{2,}')


def txt_to_hdf5(infile_name: str, outfile_name: str) -> None:
    """
    A function that reads a COMSOL .txt file and stores a sparse specially constructed HDF5 file,
    while keeping metadata from the original file.

    The constructed HDF5 file will have the following entries:
    - index: Each of these is the set of indices found for each dimension, in order of occurrence
      - x
      - y
      - z
    - data: Each of these is one full row of values, in X*Y*Z order (X and Z are flipped wrt. the original txt file!)
      - <one for each column, matching the column name from the file>

    The metadata (the lines starting with %) is stored as an attribute on the HDF5 root.
    Columns specifying a unit (separated by one space) have this information attached as metadata on data/<...>.

    :param infile_name: The input file name. Must be a COMSOL-generated .txt file, or the behaviour will be undefined.
    :param outfile_name: The output file name.
    :return: None.
    """
    f = open(infile_name, 'r')

    metadata, headers = [], []
    x, y, z = OrderedDict(), OrderedDict(), OrderedDict()
    values = {}

    prev_line = ""
    for line in f:
        if line.startswith('%'):
            metadata.append(line)
            prev_line = line
            continue
        elif prev_line.startswith('%'):
            # The previous line is the "header line" if this one does not start with a %
            headers = [tuple(x.split()) for x in re.split(gte2spaces, prev_line[2:].strip())]
            # Discard the headers for X/Y/Z
            headers = headers[3:]
            # Reset the previous line
            prev_line = ""
            # Store an empty array in values for each header name
            for header in headers:
                values[header[0]] = []

        cols = line.split()
        x[float(cols[0])] = None  # store x, y, z uniquely and in order of occurrence
        y[float(cols[1])] = None
        z[float(cols[2])] = None

        for i in range(len(headers)):
            values[headers[i][0]].append(float(cols[3+i]))

    # Construct numpy arrays from the collected X/Y/Z index set and store them
    x_index, y_index, z_index = [np.array(list(i.keys()), dtype=float) for i in [x, y, z]]

    nx, ny, nz = len(x_index), len(y_index), len(z_index)
    for i in range(len(headers)):
        val_arr = np.array(values[headers[i][0]], dtype=np.float)
        # Z iterates the slowest (see the txt files generated by COMSOL), so we reorder to get x,y,z
        # order like typical products (e.g. itertools.product) would generate. hdf5_to_dataframe will construct the
        # index in x,y,z product order, where X iterates the slowest, so we reshape the data according to this.
        val_arr = val_arr.reshape((nz, ny, nx)).swapaxes(0, 2).reshape((nz*ny*nx,))
        val_mat = csr_matrix(np.nan_to_num(val_arr))
        values[headers[i][0]] = val_mat

    # Use h5sparse to save (lots of) space when storing (lots of) NaN or 0.0 values
    with h5sparse.File(outfile_name) as h5f:
        # Store the metadata from the .txt file as attributes on the root of the hdf5 file
        h5f.attrs.metadata = "\n".join(metadata)
        # Create a dataset for each column, store the np.array and the original index of the column
        for i in range(len(headers)):
            dataset = h5f.create_dataset(f'data/{headers[i][0]}', data=values[headers[i][0]])
            dataset.attrs['column_index'] = i  # store this to preserve column order
            if len(headers[i]) > 1:
                # If the header for this column has more than one part, it's assumed the second part is the unit
                # and this is also stored as metadata
                dataset.attrs['unit'] = headers[i][1]

        # Store the index definitions for x/y/z
        h5f.create_dataset('index/x', data=x_index)
        h5f.create_dataset('index/y', data=y_index)
        h5f.create_dataset('index/z', data=z_index)


def hdf5_to_data_frame(filename: str) -> pd.DataFrame:
    """
    Reads in an HDF5 file (as written by txt_to_hdf5) and returns a pandas DataFrame matching it.
    Note that not all metadata stored in the HDF5 file is also attached to the DataFrame as of now.
    :param filename: The HDF5 file to read in. Must be written by or adhere to the format that txt_to_hdf5 defines.
    :return: A pandas DataFrame matching the data stored in the HDF5 file.
    """
    with h5sparse.File(filename, 'r') as h5f:
        x, y, z = [h5f[f'index/{i}'] for i in ['x', 'y', 'z']]
        column_indices = []
        column_headers = []
        column_values = []

        # Populate the values dict from the columns in 'data'
        for column_name in h5f['data']:
            col = h5f['data'][column_name]
            column_headers.append(column_name)
            column_indices.append(col.attrs['column_index'])
            column_values.append(col[:].toarray())

        # Reorder the column data and headers according to the (inverse) permutation gathered from the metadata
        column_indices = np.argsort(column_indices)  # argsort constructs the inverse permutation
        column_headers = np.array(column_headers, dtype=str)[column_indices]
        # Values needs to be (n_d x n_c), not (n_c x n_d) for the DataFrame constructor
        # where n_d is the number of data rows and n_c is the number of columns. So swap the axes.
        column_values = np.vstack(column_values)[column_indices].swapaxes(0, 1)

        # Generate a MultiIndex from the product of x/y/z (X iterates the 'slowest' here)
        idx = pd.MultiIndex.from_product([x, y, z], names=('x', 'y', 'z'))
        # Construct the DataFrame from the values with the index
        df = pd.DataFrame(data=column_values, columns=column_headers, index=idx)

        return df


def data_frame_to_data_grid(df: pd.DataFrame) -> Tuple[np.array, np.array, np.array, np.array]:
    """
    Turns a DataFrame (as constructed by hdf5_to_data_frame) into a data grid that can be used with the
    scipy.interpolate.RegularGridInterpolator.
    The result can be passed on like: RegularGridInterpolator((x,y,z), data_grid).
    :param df: The pandas DataFrame. Must be constructed by or adhere to the format that hdf5_to_data_frame defines.
    :return: A 4-tuple of numpy arrays: (x, y, z, data_grid).
    """
    x, y, z = map(lambda level: level.values, df.index.levels)  # expect 3 index levels or die (destructuring will fail)

    # Gather number of indices along each axis from the DataFrame's MultiIndex
    nx, ny, nz = [len(df.index.levels[i]) for i in range(len(df.index.levels))]
    # Create the data grid
    data_grid = np.nan_to_num(df.values.reshape((nx, ny, nz, -1)))

    # Return x/y/z and the data grid. Both are needed to construct a RegularGridInterpolator.
    return x, y, z, data_grid


def hdf5_to_data_grid(filename: str) -> Tuple[np.array, np.array, np.array, np.array]:
    """
    Shortcut function to read in an HDF5 file (as written by txt_to_hdf5) and turn it into a data grid.
    Refer to hdf5_to_data_frame and data_frame_to_data_grid for further info; this function is defined purely
    in terms of the composition of the two.
    :param filename: The HDF5 file to read in.
    :return: A 4-tuple of numpy arrays: (x, y, z, data_grid).
    """
    return data_frame_to_data_grid(hdf5_to_data_frame(filename))


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""