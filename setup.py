#!/usr/bin/env python3
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMInject
#
# Copyright (C) 2018,2020 CFEL Controlled Molecule Imaging group
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# If you use this program for scientific work, you should correctly reference it; see the LICENSE.md file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

import sys
from setuptools import setup, find_packages, Extension, dist

# https://luminousmen.com/post/resolve-cython-and-numpy-dependencies
dist.Distribution().fetch_build_eggs(['Cython>=0.29.15', 'numpy>=1.16.0'])
from Cython.Build import cythonize
import numpy

name = 'cminject'
version = '1.0.0'
release = version
copyright = 'Muhamed Amin, Simon Welker, and the CFEL Controlled Molecule Imaging group'

long_description = """
CMInject -- A Python framework for defining and executing particle trajectory simulations for sample injection.

Developed by Muhamed Amin, Simon Welker, and the Controlled Molecule Imaging group at the Center for
Free-Electron Laser Science, Deutsches Elektronen-Synchrotron DESY and Universität Hamburg, Hamburg,
Germany.

Original author:    Muhamed Amin <muhamed.amin@cfel.de> and the CMI COMOTION team
Current maintainer: Simon Welker <simon.welker@cfel.de> and the Controlled Molecule Imaging group
"""

if sys.version_info < (3, 6):
    sys.exit('Sorry, Python < 3.6 is not supported')

package_dir = {'': 'lib'}
packages = find_packages(where='lib')

scripts = [
    'bin/cminject',
    'bin/cminject_txt-to-hdf5',
    'bin/cminject_visualize',
    'bin/cminject_analyze-asymmetry'
]

install_requires = [
    'scipy>=1.3.0',
    'numpy>=1.16.0',
    'numba~=0.50.1',
    'pandas>=0.24.0',
    'matplotlib>=3.1.0',
    'h5py>=2.10.0',
    'sphinx>=2.4.4',
    'docutils<0.16',  # see e.g. https://github.com/matplotlib/matplotlib/pull/16358
    'sphinx_rtd_theme~=0.5.0',
    'tqdm>=4.41.1'
]

extensions = [
    Extension('cminject.utils.cython_interpolation',
              ['lib/cminject/utils/cython_interpolation.pyx'],
              include_dirs=[numpy.get_include(), '.'])
]

classifiers = [
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3 :: Only',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development :: Libraries',
]

build_sphinx_options = {
    'project': ('setup.py', name),
    'version': ('setup.py', version),
    'release': ('setup.py', release),
    'source_dir': ('setup.py', 'doc'),
    'copyright': ('setup.py', copyright)
}

setup(name=name,
      version=release,
      description='A Python framework for defining and executing particle trajectory simulations '
                  'for sample injection.',
      long_description=long_description,
      author='Simon Welker, Muhamed Amin, and the CFEL Controlled Molecule Imaging group',
      author_email='simon.welker@cfel.de',
      maintainer='CFEL Controlled Molecule Imaging group',
      maintainer_email='cminject@desy.de',  # need to set up an email list with that name ;-)
      url='https://github.com/CFEL-CMI/cminject',  # put on github and adjust
      package_dir=package_dir,
      packages=packages,
      scripts=scripts,
      python_requires='>=3.6',
      install_requires=install_requires,
      ext_modules=cythonize(extensions),
      command_options={'build_sphinx': build_sphinx_options,},
      classifiers=classifiers
)



