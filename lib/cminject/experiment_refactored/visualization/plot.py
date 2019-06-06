from typing import List

import numpy as np
from cminject.experiment_refactored.definitions.base import Particle

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_particles(particles: List[Particle], plot_trajectories=False):
    xs0, ys0, zs0 = [[p.initial_position[i] for p in particles] for i in range(3)]
    xs, ys, zs = [[p.position[i] for p in particles] for i in range(3)]
    color = [
        'red' if p.lost and not p.reached_any_detector else 'green'
        for p in particles
    ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.scatter(xs, ys, zs=zs, s=3, c=color)
    plt.scatter(xs0, ys0, zs=zs0, s=3, c='black')

    for p in particles:
        for detector_id, hits in p.detector_hits.items():
            hits = np.array([hit.hit_position for hit in hits])
            plt.scatter(hits[:, 0], hits[:, 1], zs=hits[:, 2], s=10, color='yellow')

        if plot_trajectories and p.trajectory:
            trajectory = np.array(p.trajectory)
            ts = trajectory[:, 0]
            ps = trajectory[:, 1:4]
            xs = ps[:, 0]
            ys = ps[:, 1]
            zs = ps[:, 2]

            vs = trajectory[:, 4:]
            plt.plot(xs, ys, zs=zs, color='grey')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()