#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 100 -*-
#
# Copyright (C) 2008,2017 Jochen Küpper <jochen.kuepper@cfel.de>


import os
from setuptools import setup

extra_compile_args = []
library_dirs = []

long_description = """CMI RichMol -- time-dependent wave packet propagation


Python extensions for Stark effect calculations of molecules.

Developed by Andrey Yachmenev and the Controlled Molecule Imaging Group at the Center for
Free-Electron Laser Science, Deutsches Elektronen-Synchrotron DESY and Universität Hamburg, Hamburg,
Germany.

Original author:    Andrey Yachmenev <andrey.yachmenev@cfel.de>
Current maintainer: Andrey Yachmenev <andrey.yachmenev@cfel.de>
See the distribution files AUTHORS and THANKS for further contributions.
"""


setup(name="richmol",
      author              = "Andrey Yachmenev and the CFEL-CMI group",
      author_email        = "andrey.yachmenev@cfel.de",
      maintainer          = "Andrey Yachmenev and the CFEL-CMI group",
      maintainer_email    = "andrey.yachmenev@cfel.de",
      url                 = "https://controlled-molecule-imaging.org/software/richmol",
      description         = "CMI richmol",
      version             = "0.0.dev1",
      long_description    = long_description,
      license             = "GPL",
      packages            = ['lib'],
      scripts             = ['bin/richmol'],
      test_suite          = 'tests',
      )
