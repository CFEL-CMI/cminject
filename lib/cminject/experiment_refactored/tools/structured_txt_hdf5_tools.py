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
from typing import Tuple, List

import h5py
import h5sparse

from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd

import re

from collections import OrderedDict

gte2spaces = re.compile(r'\s{2,}')

index_names = ['x', 'y', 'z']


def txt_to_hdf5(infile_name: str, outfile_name: str, dimensions: int = 3) -> None:
    """
    A function that reads a structured .txt file and stores a sparse specially constructed HDF5 file,
    while keeping metadata from the original file. The general structure of the .txt input file must be as follows:

    ``
    % <Some meta information>
    % <Some more meta information>
    % X  Y  Z  A (m)  B (s)  C ...
    0.0  0.0  0.0  1.0  2.0  3.0 ...
    0.1  0.0  0.0  4.0  5.0  6.0 ...
    ...
    ``

    This matches the .txt output format of COMSOL for 3D grids.

    All but the last line starting with % at the beginning of the file are considered to be a general textual
    description of the file (metadata about the file itself), and are stored.

    The last line starting with % is considered the header of the columns, where each column name has to be separated
    by at least two spaces from the other column names. A single space is considered to still be part of the same
    column, and the second part of the column name is - if present - considered to be the unit of the column.
    This will be stored as metadata too.

    The number of columns after X/Y/Z is arbitrary.
    NOTE: X is the dimension that increments first. This convention MUST be adhered to.

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
    :param dimensions: The number of spatial dimensions in the txt file. 3 by default, will be 2 or 3 for most cases.
    :return: None.
    """
    f = open(infile_name, 'r')
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
        print("WARNING: The .txt file to convert has no headers, and so we can't assign names for columns. "
              "If you know the relevant info, it's best to add these headers by hand and run the conversion again.")
        metadata = []
        # Read the first line and rewind again
        first_line = f.readline()
        f.seek(0)
        # Get the number of columns by splitting the line (on whitespace), and construct the header as just 0-based ints
        headers = [(i,) for i in range(len(first_line.split()) - dimensions)]
        for header in headers:
            values[header[0]] = []

    # From here on we assume every line to only contain data
    for line in f:
        # After the header handling: Split current line on spaces, update index and values
        cols = line.split()
        for d in range(dimensions):
            # Add an entry to the index for the current row, for each dimension
            index[d][float(cols[d])] = None  # We only use the key, value can be None
        for i in range(len(headers)):
            values[headers[i][0]].append(float(cols[dimensions+i]))

    # Construct numpy arrays for each collected index set per dimension
    index = [np.array(list(i.keys()), dtype=float) for i in index]
    index_n = [len(i) for i in index]

    for i in range(len(headers)):
        val_arr = np.array(values[headers[i][0]], dtype=np.float)
        # Transpose the whole array to match "normal" X/Y/Z/... iteration order. The order that the dimensions
        # increment in is inverted in the txt file format, so transpose it here.
        val_arr = val_arr.reshape(tuple(reversed(index_n))).transpose().reshape((np.prod(index_n),))
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

        # Store the index definitions for each spatial dimension
        for i in range(dimensions):
            h5f.create_dataset(f'index/{i}', data=index[i])


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
            column_values.append(col[:].toarray())

        # Reorder the column data and headers according to the (inverse) permutation gathered from the metadata
        column_indices = np.argsort(column_indices)  # argsort constructs the inverse permutation
        column_headers = np.array(column_headers, dtype=str)[column_indices]
        # Values needs to be (n_d x n_c), not (n_c x n_d) for the DataFrame constructor
        # where n_d is the number of data rows and n_c is the number of columns. So swap the axes.
        column_values = np.vstack(column_values)[column_indices].swapaxes(0, 1)

        # Generate a MultiIndex from the product of the indices for each spatial dimension
        idx = pd.MultiIndex.from_product(index, names=tuple(index_names[:dimensions]))
        # Construct the DataFrame from the values with the index
        df = pd.DataFrame(data=column_values, columns=column_headers, index=idx)

        return df


def data_frame_to_data_grid(df: pd.DataFrame) -> List[np.array]:
    """
    Turns a DataFrame (as constructed by hdf5_to_data_frame) into a data grid that can be used with the
    scipy.interpolate.RegularGridInterpolator.
    The result can be passed on like: RegularGridInterpolator((x,y,z), data_grid).
    :param df: The pandas DataFrame. Must be constructed by or adhere to the format that hdf5_to_data_frame defines.
    :return: A (n+1)-tuple of numpy arrays: (x, [y, [z, ...]], data_grid), where n is the number of spatial dimensions
    """
    # Gather index and number of indices along each axis from the DataFrame's MultiIndex
    index = [level.values for level in df.index.levels]
    index_n = [len(df.index.levels[i]) for i in range(len(df.index.levels))]

    # Create the data grid and shape it into expected n-dimensional form
    new_shape = tuple(index_n + [-1])
    data_grid = np.nan_to_num(df.to_numpy().reshape(new_shape))
    # Return the index array and the data grid. Both are needed to construct a RegularGridInterpolator.
    return index + [data_grid]


def hdf5_to_data_grid(filename: str) -> List[np.array]:
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