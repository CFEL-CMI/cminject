import os
import matplotlib.pyplot as plt

from cminject.result_storages import HDF5ResultStorage
from cminject.utils.visualization import plot_detectors, plot_trajectories_colored

with HDF5ResultStorage('example_output.h5', 'r') as f:
    plot_trajectories_colored(f.get_trajectories())
    plt.savefig('example_trajectories.pdf')

with HDF5ResultStorage('als_output.h5', 'r') as f:
    detectors = f.get_detectors()
    plot_detectors(detectors, 'x,vx', bins=(20, 20))
    plt.savefig('als_detectors_x_vx.pdf')
    plot_detectors(detectors, 'vx,vz', bins=(20, 20))
    plt.savefig('als_detectors_vx_vz.pdf')

plt.show()
os.system("cminject_analyze-beam als_output.h5 -d als_experimental_focus_curve.txt -o als_experiment_comparison.pdf")
