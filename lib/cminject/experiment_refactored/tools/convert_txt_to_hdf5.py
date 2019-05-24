#!/usr/bin/env python3
"""
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
import argparse

from cminject.experiment_refactored.tools.comsol_hdf5_tools import txt_to_hdf5


def main():
    parser = argparse.ArgumentParser(prog='convert_txt_to_hdf5',
                                     formatter_class=argparse.MetavarTypeHelpFormatter)
    parser.add_argument('-i', help='Input file (COMSOL .txt format)', type=str, required=True)
    parser.add_argument('-o', help='Output file (.hdf5)', type=str, required=True)
    args = parser.parse_args()

    txt_to_hdf5(args.i, args.o)
    print("Done.")


if __name__ == '__main__':
    main()


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""