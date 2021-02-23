#!/bin/bash
export PYTHONPATH=".:$PYTHONPATH"
`# Use (the default OneFlowField setup with) als_field.h5`
cminject -f "als_field.h5"\
    `# Gold particles: Use appropriate density`\
    -rho 19320\
    `# Use 2 spatial dimensions`\
    -D 2\
    `# Normal distribution for x, constant for z position`\
    -p G[0,0.002] -0.128\
    `# Normal distributions for velocities in x and z`\
    -v G[0,1] G[1,0.1]\
    `# Normal distribution for the particle radius`\
    -r G[1.35e-8,1.91e-9]\
    `# Simulate the timespan from 0 to 1 seconds`\
    -t 0 1\
    `# Use time step of 10us`\
    -ts 1e-5\
    `# Insert the detectors at positions as in the experiment`\
    -d 0.001 0.00141 0.00181 0.00221 0.00261\
       0.00301 0.00341 0.00381 0.00421 0.00461\
    `# Simulate 10^3 particles`\
    -n 1000\
    `# Write results to als_output.h5`\
    -o als_output.h5\
    `# Turn Brownian motion on for this simulation`\
    -B
