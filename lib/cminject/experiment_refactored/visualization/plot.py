from typing import List

import numpy as np
from cminject.experiment_refactored.definitions.base import Particle

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_2d_from_hdf5(hd5file: str, dim_a: int, dim_b: int) -> None:
    import h5py
    with h5py.File(hd5file) as h5f:
        particles = h5f['particles']
        hits = h5f['detector_hits']

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(particles[:, dim_a], particles[:, dim_b])

        for particle_id in hits.keys():
            ax.scatter(hits[particle_id][:, dim_a], hits[particle_id][:, dim_b])

    plt.show()


if __name__ == '__main__':
    import sys
    plot_2d_from_hdf5(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))


def plot_particles(particles: List[Particle], plot_trajectories=False, dimensions: int = 3):
    if dimensions not in [2, 3]:
        raise ValueError("Can only plot for 2 or 3 dimensions.")

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
        plt.scatter(*positions, s=3, c=color)
        plt.scatter(*initial_positions, s=3, c='black')

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
    plt.show()
