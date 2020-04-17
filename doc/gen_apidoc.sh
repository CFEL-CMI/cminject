#!/usr/bin/env bash

# Generate Apidoc from cminject/ subfolder, output to ..docs/_modules/,
# ignore LBM and alignment subdirectories, do not generate toctree file
pushd ../lib/ >/dev/null
sphinx-apidoc -T -o ../doc/_modules/ cminject/ cminject/LBM cminject/alignment
popd >/dev/null