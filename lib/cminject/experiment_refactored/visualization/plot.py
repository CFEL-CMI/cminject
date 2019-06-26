from typing import List, Tuple

import h5py
import numpy as np
from cminject.experiment_refactored.definitions.base import Particle

from matplotlib import pyplot as plt


def _find_dims(description, dims):
    return [np.where(description == dim)[0][0] for dim in dims]


def get_hist2d_figure_from_hdf5(hdf5_filename: str, dimension_pairs: List[Tuple[str, str]]):
    with h5py.File(hdf5_filename) as h5f:
        detectors = h5f['detector_hits']
        detector_count = len(detectors)
        detector_keys = list(detectors.keys())
        plot_count = len(dimension_pairs)

        fig, axes = plt.subplots(detector_count, plot_count)
        for j in range(plot_count):
            for i in range(detector_count):
                hits = detectors[detector_keys[i]]

                dims = dimension_pairs[j]
                dim_a, dim_b = _find_dims(hits.attrs['description'], dims)
                ax = axes[i, j]

                x, y = hits[:, dim_a], hits[:, dim_b]
                ax.hist2d(x, y, bins=100, cmap='viridis')
                ax.autoscale(enable=True, tight=False)
                ax.set_xlabel(dims[0])
                ax.set_ylabel(dims[1])

                ax.autoscale(enable=False)
                ax.scatter([0.0], [0.0], s=100, c='r', marker='o')

        return fig, axes


def get_traj_figure_from_particles(particles: List[Particle], plot_trajectories=False):
    from mpl_toolkits.mplot3d import Axes3D
    dimensions = int(len(particles[0].position) / 2)
    if dimensions not in [2, 3]:
        raise ValueError("Can only plot 2D/3D!")

    initial_positions = np.array([p.initial_position[:dimensions] for p in particles]).transpose()
    positions = np.array([p.position[:dimensions] for p in particles]).transpose()
    color = [
        'red' if p.lost and not p.reached_any_detector else 'green'
        for p in particles
    ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d') if dimensions == 3 else fig.add_subplot(111)

    if dimensions == 3:
        plt.scatter(positions[0], positions[1], zs=positions[2], s=3, c=color)
        plt.scatter(initial_positions[0], initial_positions[1], zs=initial_positions[2], s=3, c='black')
    elif dimensions == 2:
        plt.scatter(positions[0], positions[1], s=3, c=color)
        plt.scatter(initial_positions[0], initial_positions[1], s=3, c='black')

    for p in particles:
        for detector_id, hits in p.detector_hits.items():
            hits = np.array([hit.hit_position for hit in hits])
            if dimensions == 3:
                plt.scatter(hits[:, 0], hits[:, 1], zs=hits[:, 2], s=10, color='yellow')
            elif dimensions == 2:
                plt.scatter(hits[:, 0], hits[:, 1], s=10, color='yellow')

        if plot_trajectories and p.trajectory:
            trajectory = np.array(p.trajectory)
            ts = trajectory[:, 0]
            positions = trajectory[:, 1:(dimensions+1)].transpose()
            plt.plot(*positions, color='grey')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if dimensions == 3:
        ax.set_zlabel('z')

    return fig