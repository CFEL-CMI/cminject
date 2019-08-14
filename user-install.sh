#!/usr/bin/env bash
python3 -m pip install --user Cython numpy  # install requirements for setup.py
python3 setup.py build_ext --inplace        # use setup.py to build extensions inplace
python3 setup.py install --user             # install the software, which also installs dependencies
