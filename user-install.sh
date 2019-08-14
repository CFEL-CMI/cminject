#!/usr/bin/env bash

# Try importing cython before installing it
python3 -c "import Cython; import cython"
if [[ $? -ne 0 ]]; then
    echo "Cython not found, installing to user dir..." >&2
    python3 -m pip install --user Cython
fi

# Try importing numpy before installing it
python3 -c "import numpy"
if [[ $? -ne 0 ]]; then
    echo "numpy not found, installing to user dir..." >&2
    python3 -m pip install --user numpy
fi

python3 setup.py build_ext --inplace        # use setup.py to build extensions inplace
python3 setup.py install --user             # install the software, which also installs dependencies
