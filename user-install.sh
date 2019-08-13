#!/usr/bin/env bash

# install CMI injector in current user's home directory
python3 setup.py build_ext --inplace
python3 setup.py install --user
