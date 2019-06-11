from typing import List

import numpy as np
from cminject.experiment_refactored.definitions.base import Particle

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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
