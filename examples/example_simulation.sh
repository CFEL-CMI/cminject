#!/bin/bash
export PYTHONPATH=".:$PYTHONPATH"
cminject \
  `# Use ExampleSetup from example_setup.py` \
  -s example_setup.ExampleSetup \
  `# Use the 2D example field` \
  -f example_field.h5 \
  `# Simulate 100 particles` \
  -n 100 \
  `# Write the results to example-output.h5` \
  -o example_output.h5 \
  `# Track and store trajectories` \
  -T
