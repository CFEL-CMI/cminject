#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np
from cminject.result_storages import HDF5ResultStorage


# Global settings for LaTex output
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = '\DeclareUnicodeCharacter{2212}{-}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 16


# Create figure and axes
fig = plt.figure(figsize=(6, 4))
ax = fig.gca()


# Simulation
filename = sys.argv[1] if len(sys.argv) > 1 else './als_output.h5'
with HDF5ResultStorage(filename) as rs:
    detectors = rs.get_detectors()

    # Get x positions and retrieve 70th percentiles of x distribution for each detector
    hits = np.array(list(detectors.values()))
    xs = hits['position'][:, :, 0]
    x_percentiles = 2 * np.percentile(np.abs(xs), 70, axis=1)

    # Get z positions as the last dimension of the first entry of each detector
    # This is reasonable to do since we know our detectors detect at exactly one Z position
    zs = hits['position'][:, 0, -1]

    # Plot accordingly scaled values (z in mm, x percentiles in um)
    ax.plot(zs * 1e3, x_percentiles * 1e6, '-o', label='Simulation', zorder=-1)
    xs_ = np.abs(xs)


# Experiment
x_spheres = [1.41, 1.81, 2.21, 2.61, 3.01, 3.41, 3.81, 4.21, 4.61]
y_spheres = [49.9, 42.6, 38.3, 35.8, 35.8, 44.44, 50.62, 64.2, 69.1]
yerr_spheres = [1.84, 1.85, 1.74, 3.49, 3.49, 3.03, 3.49, 6.98, 9.24]
ax.errorbar(x_spheres, y_spheres, yerr_spheres, label='Experiment')


# Labeling, x/y limits and grid
ax.legend(loc='lower right')
ax.set_xlabel('$z$ distance from ALS exit [mm]')
ax.set_ylabel('70th percentile of $x$ positions [$\mu$m]')
ax.set_xlim(1.1, 4.9)
ax.set_ylim(0, 100)
ax.grid(axis='y')


# Show and save figure
plt.show()
fig.savefig('focus-curve-gold-27nm-uniform-zoom.pdf', bbox_inches="tight")
