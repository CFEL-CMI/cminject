#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 100 -*-
#
# Copyright (C) 2008,2017 Jochen Küpper <jochen.kuepper@cfel.de>


import os
from setuptools import setup
from setuptools import setup, find_packages

long_description = """CMI Injector -- Simulating particles' trajectories in different forcefields


Developed by Muhamed Amin and the Controlled Molecule Imaging Group at the Center for
Free-Electron Laser Science, Deutsches Elektronen-Synchrotron DESY and Universität Hamburg, Hamburg,
Germany.

Original author:    Muhamed Amin <muhamed.amin@cfel.de>
Current maintainer: Muhamed Amin <muhamed.amin@cfel.de>
"""

package_dir = {"": "lib"}

packages = find_packages(where="lib")

provides = [
    'cmiinject',
]


requires = [
    'Python (>=3.0)',
    'scipy (>=0.14.1)',
    'numpy (>= 1.7.1)']


setup(name="cmi-inject",
      author              = "Muhamed Amin and the CFEL-CMI group",
      author_email        = "muhamed.amin@cfel.de",
      maintainer          = "muhamed amin and the CFEL-CMI group",
      maintainer_email    = "muhamed.amin@cfel.de",
      url                 = "https://stash.desy.de/projects/CMIFLY/repos/cmi-injector/browse",
      description         = "CMI inject",
      version             = "0.0.0",
      long_description    = long_description,
      license             = "GPL",
      package_dir         = package_dir,
      packages            = packages,
      scripts             = None,
      requires            = requires,
      install_requires    = requires,
      )

