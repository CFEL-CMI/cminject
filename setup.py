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
from setuptools import setup, find_packages, Extension

# ---- Project metadata ----
name = 'cminject'
version = '1.0.8'
release = version
author = 'Simon Welker, Surya Kiran Peravali, Muhamed Amin, and the CFEL Controlled Molecule Imaging group'
curr_maintainer_name = 'Surya Kiran Peravali'
curr_maintainer_email = 'surya.kiran.peravali@desy.de'
description = 'A Python framework for defining and executing particle trajectory simulations for sample injection.'

long_description = f"""
CMInject -- {description}.

Developed by {author} at the Center for Free-Electron Laser Science, Deutsches Elektronen-Synchrotron DESY and
Universität Hamburg, Hamburg, Germany.

Original author:    Muhamed Amin <muhamed.amin@cfel.de> and the CMI COMOTION team
Current maintainer: {curr_maintainer_name} <{curr_maintainer_email}> and the Controlled Molecule Imaging group
"""

if sys.version_info < (3, 6):
    sys.exit('Sorry, Python < 3.6 is not supported')

# ---- Package discovery ----
package_dir = {'': 'lib'}
packages = find_packages(where='lib')

scripts = [
    'bin/cminject',
    'bin/cminject_txt-to-hdf5',
    'bin/cminject_visualize',
    'bin/cminject_analyze-asymmetry',
    'bin/cminject_analyze-beam'
]

# ---- Runtime dependencies (version-conditioned) ----
# Keep legacy pins for Py<3.9, use modern wheels for Py>=3.9
install_requires = [
    # numpy
    'numpy>=1.23.5,<2.0; python_version >= "3.9"',
    'numpy>=1.19.5,<1.20; python_version < "3.9"',

    # scipy
    'scipy>=1.10,<1.12; python_version >= "3.9"',
    'scipy>=1.5,<1.6; python_version < "3.9"',

    # numba
    'numba>=0.57,<0.59; python_version >= "3.9"',
    'numba>=0.53,<0.54; python_version < "3.9"',

    # pandas
    'pandas>=1.5.3,<2.0; python_version >= "3.9"',
    'pandas>=1.0.5,<1.1; python_version < "3.9"',

    # matplotlib
    'matplotlib>=3.7,<3.9; python_version >= "3.9"',
    'matplotlib>=3.3.4,<3.4; python_version < "3.9"',

    # h5py
    'h5py>=3.8,<3.11; python_version >= "3.9"',
    'h5py>=3.1,<3.2; python_version < "3.9"',

    # tqdm
    'tqdm>=4.66,<5',
]

# ---- Optional deps (docs moved out of install_requires) ----
extras_require = {
    'docs': [
        # modern docs toolchain for Py>=3.9
        'Sphinx>=7,<8; python_version >= "3.9"',
        'sphinx_rtd_theme>=2.0,<3; python_version >= "3.9"',
        'docutils>=0.19,<0.21; python_version >= "3.9"',

        # legacy docs toolchain for Py<3.9
        'Sphinx~=3.4.3; python_version < "3.9"',
        'sphinx_rtd_theme~=0.5.1; python_version < "3.9"',
        'docutils<0.16; python_version < "3.9"',
    ],
    'dev': [
        'pytest>=7,<8',
    ],
}

# ---- Extensions ----
# Note: NumPy headers are provided by pyproject.toml (build-system.requires)
try:
    import numpy
except Exception:
    # setuptools with build isolation will install numpy before reaching here
    numpy = None

include_dirs = []
if numpy is not None:
    include_dirs.append(numpy.get_include())
include_dirs.append('.')

extensions = [
    Extension(
        'cminject.utils.cython_interpolation',
        ['lib/cminject/utils/cython_interpolation.pyx'],
        include_dirs=include_dirs,
        language='c'
    )
]

# cythonize at build time (requires cython in pyproject)
from Cython.Build import cythonize
ext_modules = cythonize(
    extensions,
    language_level=3,
    compiler_directives={"boundscheck": False, "wraparound": False}
)

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
    'Programming Language :: Python',
    'Programming Language :: Python :: 3 :: Only',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development :: Libraries',
]

build_sphinx_options = {
    'project': ('setup.py', name),
    'version': ('setup.py', version),
    'release': ('setup.py', release),
    'source_dir': ('setup.py', 'doc'),
    'copyright': ('setup.py', f'2021, {author}')
}

setup(
    name=name,
    version=release,
    description=description,
    long_description=long_description,
    author=author,
    author_email=curr_maintainer_email,
    maintainer='CFEL Controlled Molecule Imaging group',
    maintainer_email='cminject@desy.de',
    url='https://github.com/CFEL-CMI/cminject',
    package_dir=package_dir,
    packages=packages,
    scripts=scripts,
    python_requires='>=3.6',
    install_requires=install_requires,
    extras_require=extras_require,
    # REMOVE legacy setup_requires/fetch_build_eggs (deprecated)
    tests_require=['pytest~=6.2.3'],
    ext_modules=ext_modules,
    command_options={'build_sphinx': build_sphinx_options},
    classifiers=classifiers
)

