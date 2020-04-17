#!/usr/bin/env python3
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

import sys

from setuptools import setup, find_packages, Extension
try:
    from Cython.Build import cythonize
except ImportError:
    print("FATAL: CMInject setup requires Cython to be installed, since it has Cython extensions defined.\n"
          "Please install the 'Cython' package in your Python environment to proceed.", file=sys.stderr)
    sys.exit(1)
try:
    import numpy
except ImportError:
    print("FATAL: CMInject setup requires numpy to be installed, since its Cython extensions depend on it.\n"
          "Please install the 'numpy' package in your Python environment to proceed.", file=sys.stderr)
    sys.exit(1)


long_description = """CMI Injector -- Simulating particles' trajectories in different forcefields

Developed by Muhamed Amin, Simon Welker and the Controlled Molecule Imaging Group at the Center for
Free-Electron Laser Science, Deutsches Elektronen-Synchrotron DESY and Universität Hamburg, Hamburg,
Germany.

Original author:    Muhamed Amin <muhamed.amin@cfel.de>
Current maintainer: Muhamed Amin <muhamed.amin@cfel.de>
"""

if sys.version_info < (3, 6):
    sys.exit('Sorry, Python < 3.6 is not supported')

package_dir = {'': 'lib'}
packages = find_packages(where='lib')

scripts = [
    'bin/cminject',
    'bin/cminject_txt-to-hdf5',
    'bin/cminject_visualize',
    'bin/cminject_reconstruct-detectors',
    'bin/cminject_analyze-asymmetry',
    'bin/cminject_analyze-beam'
]

install_requires = [
    'scipy~=1.3',
    'numpy~=1.16',
    'numba~=0.44',
    'pandas~=0.24',
    'matplotlib~=3.1',
    'h5py~=2.9',
    'h5sparse~=0.1',
    'sphinx~=2.1',
    'sphinx_rtd_theme~=0.4',
    'tqdm~=4.41',
    'scikit-learn~=0.22'
]

extensions = [
    Extension('cminject.utils.interpolation.cython_interpolation',
              ['lib/cminject/utils/interpolation/cython_interpolation.pyx'],
              include_dirs=[numpy.get_include(), '.'])
]

name = 'cminject'
version = '0.1'
release = '0.1.0'

build_sphinx_options = {
    'project': ('setup.py', name),
    'version': ('setup.py', version),
    'release': ('setup.py', release),
    'source_dir': ('setup.py', 'doc')
}

setup(
    name=name,
    version=release,
    description='A framework for particle injection trajectory simulations in different force fields',
    long_description=long_description,
    author='Simon Welker, Muhamed Amin, and the CFEL-CMI group',
    author_email='simon.welker@cfel.de',
    maintainer='CFEL-CMI group',
    maintainer_email='cminject@cfel.de',  # need to set up an email list with that name ;-)
    url='https://stash.desy.de/projects/CMIFLY/repos/cmi-injector/browse',
    package_dir=package_dir,
    packages=packages,
    scripts=scripts,
    python_requires='>=3.7',
    install_requires=install_requires,
    setup_requires=['Cython>=0.29.10', 'numpy>=1.16.0'],
    ext_modules=cythonize(extensions),
    command_options={
        'build_sphinx': build_sphinx_options,
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
    ],
)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
