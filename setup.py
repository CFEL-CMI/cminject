#!/usr/bin/env python3.4
# -*- coding: utf-8; fill-column: 100 -*-
#
# Copyright (C) 2008,2017 Jochen Küpper <jochen.kuepper@cfel.de>


import sys

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

from setuptools import find_packages

long_description = """CMI Injector -- Simulating particles' trajectories in different forcefields

Developed by Muhamed Amin, Simon Welker and the Controlled Molecule Imaging Group at the Center for
Free-Electron Laser Science, Deutsches Elektronen-Synchrotron DESY and Universität Hamburg, Hamburg,
Germany.

Original author:    Muhamed Amin <muhamed.amin@cfel.de>
Current maintainer: Muhamed Amin <muhamed.amin@cfel.de>
"""

if sys.version_info < (3,6):
    sys.exit('Sorry, Python < 3.6 is not supported')


package_dir = {"": "lib"}
packages = find_packages(where="lib")

provides = [
    'cminject',
    'txt_to_hdf5',
    'visualize'
]

install_requires = [
    'scipy>=1.3.0',
    'numpy>=1.16.0',
    'pandas>=0.24.0',
    'matplotlib>=3.1.0',
    'h5py>=2.9.0',
    'h5sparse>=0.1.0',
    'sphinx_rtd_theme~=0.4.3'
]

extensions = [
    Extension("cminject.utils.cython_interpolation",
              ["lib/cminject/utils/cython_interpolation.pyx"],
              include_dirs=[numpy.get_include(), "."])
]

setup(
    name="cminject",
    author="Simon Welker, Muhamed Amin and the CFEL-CMI group",
    author_email="simon.welker@cfel.de",
    maintainer="Muhamed Amin and the CFEL-CMI group",
    maintainer_email="muhamed.amin@cfel.de",
    url="https://stash.desy.de/projects/CMIFLY/repos/cmi-injector/browse",
    description="A framework for particle injection trajectory simulations in different force fields",
    version="0.1.0",
    long_description=long_description,
    license="GPL",
    package_dir=package_dir,
    packages=packages,
    scripts=None,
    provides=provides,
    python_requires='>=3.7',
    install_requires=install_requires,
    ext_modules=cythonize(extensions)
)
