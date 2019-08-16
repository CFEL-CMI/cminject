#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of CMInject
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this program for scientific work, you should correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

from abc import abstractmethod, ABC
from typing import List, Tuple, Any, Callable, Union

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter


class Visualizer(ABC):
    """
    A generic class that should be implemented to do different kinds of visualization of results.

    *Always* based on a filename and not on in-memory objects, due to the philosophy that everything you'll visualize
    is part of your result, so it should be stored as data in your result file, and you should be able to construct
    a visualization based solely on that result data. This is also more sensible for sharing data and visualizations
    with other people and looking at previously done visualizations again.
    """
    def __init__(self, filename: str):
        self.filename = filename

    @abstractmethod
    def visualize(self) -> Any:
        """
        A highly generic method to construct a visualization based on the filename passed at construction.

        Note that it's not conceptually required that calling this method shows anything on screen, see e.g.
        the docstring for MatplotlibVisualizer as to why this is.

        :return: Any set of objects ne. See e.g. MatplotlibVisualizer for a more concretely types
        """
        pass


class MatplotlibVisualizer(Visualizer):
    """
    A more concrete visualizer class (with regards to the expected return type of `visualize`).

    Subclasses are expected to return a Matplotlib figure and a set of axes, which the caller is able to do
    additional modifications (labeling, sizing, ...) on before actually showing the visualization.
    """
    @abstractmethod
    def visualize(self) -> Tuple[Figure, np.ndarray]:
        pass


class DetectorHistogramVisualizer(MatplotlibVisualizer):
    """
    A class to visualize 1D/2D histograms on a detector for arbitrary dimension pairs,
    based on an HDF5 file and Matplotlib.

    Whether a 1D or 2D histogram is displayed for each dimension pair is determined based on the data in both
    dimensions: If either of the value arrays holds the same value for every point, it's not used and a 1D histogram
    of the other dimension is drawn in its place. In all other cases, 2D histograms are drawn. (If the value arrays
    are the same along both dimensions, a 1D histogram for the second dimension will be drawn, but it will be rather
    boring).
    """
    # describes a dimension as either a string or a callable that takes a list of np.arrays and returns one np.array
    # (a function that works on a list of hits and returns an array of values, one per hit)
    DimensionType = Union[str, Callable[[List[np.array]], np.array]]

    def __init__(self, filename: str, dimension_pairs: List[Tuple[DimensionType, DimensionType]],
                 colormap: str = 'viridis', highlight_origin: bool = True,
                 bins: int = 100, bins_1d: int = 30):
        """
        The constructor for DetectorHistogramVisualizer.

        :param filename: The HDF5 file to read in. The expected format is the one like written by HDF5ResultStorage.
        :param dimension_pairs: A List of tuples, each defining a pair of dimensions that should be plotted against
        each other on a 2D histogram. Example: [('x', 'y'), ('x', 'z'), ('y', 'z')]. The dimensions must be present
        and annotated with their names in the HDF5 file.
        :param colormap: The matplotlib colormap to use for the histograms. 'viridis' by default.
        :param highlight_origin: A flag that determines whether (0.0, 0.0) will be 'highlighted', i.e. plotted as a
        red dot in all histograms. Useful for comparing plots along one column. True by default.
        :param bins: The number of bins in the 2d histograms, passed along to plt.hist2d().
        :param bins_1d: The number of bins in the 1d histograms, passed along to plt.hist().
        """
        super().__init__(filename)
        self.dimension_pairs = dimension_pairs
        self.colormap = colormap
        self.highlight_origin = highlight_origin
        self.bins = bins
        self.bins_1d = bins_1d

    @staticmethod
    def _get_dims(description: List[str],
                  dims: List[DimensionType],
                  hits: List[np.array])\
            -> List[np.array]:
        return [
            dim(hits) if callable(dim)
            else hits[np.where(description == dim)[0][0]]
            for dim in dims
        ]

    def _vis1d(self, ax, y, label) -> None:
        ax.hist(y, bins=self.bins_1d)
        ax.set_xlabel(label)
        ax.set_ylabel('# occurrences')

        mu, median, sigma = y.mean(), np.median(y), y.std()
        textstr = '\n'.join((r'$\mu=%.2e$' % (mu,),
                             r'$\mathrm{\overline{' + label + '}}=%.2E$' % (median,),
                             r'$\sigma=%.2E$' % (sigma,)))

        xfmt = ScalarFormatter()
        xfmt.set_powerlimits((-2, 2))
        ax.xaxis.set_major_formatter(xfmt)
        props = dict(facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, verticalalignment='top', bbox=props)

    def _vis2d(self, ax, x, y, xlabel, ylabel, cmap='viridis'):
        ax.hist2d(x, y, bins=self.bins, cmap=cmap)
        ax.set_xlabel(xlabel)  # set the labels for the 'x'/'y' axis as passed in the dimension pair
        ax.set_ylabel(ylabel)

        mu_x, median_x, sigma_x = x.mean(), np.median(x), x.std()
        mu_y, median_y, sigma_y = y.mean(), np.median(y), y.std()
        textstr = '\n'.join((
            (r'$\mu_{' + xlabel + r'}=%.2e, \mu_{' + ylabel + r'}=%.2e$') % (mu_x, mu_y),
            (r'$\mathrm{\overline{' + xlabel + r'}}=%.2E, \mathrm{\overline{' + ylabel + r'}}=%.2E$') % (median_x, median_y),
            (r'$\sigma_{' + xlabel + '}=%.2E, \sigma_{' + ylabel + '}=%.2E$') % (sigma_x, sigma_y)
        ))

        fmt = ScalarFormatter()
        fmt.set_powerlimits((-2, 2))
        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)
        props = dict(facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, verticalalignment='top', bbox=props)

    def visualize(self) -> Tuple[Figure, np.ndarray]:
        """
        Visualize the results as a plot of (potentially multiple) 1D/2D histograms.
        There will be m*n subplots arranged on a grid, with m being the number of detectors in the file,
        and n being the number of dimension pairs given at construction of this instance.

        :return: A Matplotlib figure and a 2-dimensional ndarray of all axes present in the figure.
        """
        with h5py.File(self.filename, 'r') as h5f:
            detectors = h5f['detector_hits']
            detector_count = len(detectors)
            detector_keys = list(detectors.keys())
            plot_count = len(self.dimension_pairs)

            fig, axes = plt.subplots(detector_count, plot_count)
            axes = axes.reshape((detector_count, plot_count))
            # Plot a hist2d in each subplot in the grid :)
            for j in range(plot_count):
                for i in range(detector_count):
                    # Get all hits for this detector
                    hits = detectors[detector_keys[i]]
                    description = hits.attrs['description']
                    hits = hits[:].transpose()

                    # Find the dimensions' indices and get the correct axis to plot on
                    dim_a, dim_b = self.dimension_pairs[j]
                    try:
                        x, y = self._get_dims(description, [dim_a, dim_b], hits=hits)
                    except IndexError as e:
                        raise ValueError(f"Could not find dimension {dim_a} or {dim_b}!")

                    ax = axes[i, j]
                    ax.autoscale(enable=True, tight=False)  # enable auto-scaling here :)

                    one_dimensional = False
                    if (x == x[0]).all():  # X values are all equal, plot 1D histogram of Y
                        one_dimensional = True
                        self._vis1d(ax, y, dim_b)
                    elif (y == y[0]).all():  # Y values are all equal, plot 1D histogram of X
                        one_dimensional = True
                        self._vis1d(ax, x, dim_a)
                    else:
                        self._vis2d(ax, x, y, dim_a, dim_b)

                    if self.highlight_origin and not one_dimensional:
                        # Disable auto-scaling here: If 0.0,0.0 is outside the existing bounds, that's fine, we don't
                        # want to stretch the existing plot just to show the origin marker
                        ax.autoscale(enable=False)
                        ax.scatter([0.0], [0.0], s=10, c=[(1.0, 0.0, 0.0, 0.5)], marker='x')

            return fig, axes


class TrajectoryVisualizer(MatplotlibVisualizer):
    """
    Visualizes the initial and final positions of all particles stored in an HDF5 file, along with their trajectories,
    in a 2D or 3D plot (based on the dimensionality of the stored particle positions). Detector hits are also shown.

    Useful to get a general idea of where the particles started and went, but not very useful for data analysis:
    the plots can get rather chaotic and slow (especially in 3D).
    """
    def _get_trajectories(self, h5f: h5py.File, return_velocities: bool, dimensions: int)\
            -> Union[Tuple[np.array, np.array], np.array]:
        """
        Gets the trajectory positions, and optionally velocities, from a file given a space dimensionality
        :param h5f: The h5py.File instance to retrieve the data from
        :param return_velocities: Whether to also return the velocities or only the positions
        :param dimensions: The number of spatial dimensions
        :return: Either one np.array (positions) or a 2-tuple of np.arrays (positions, velocities)
        """
        trajectories = []
        velocities = []

        for particle_id in h5f['particles']:
            if 'trajectory' in h5f[f'particles/{particle_id}'].keys():
                trajectory = h5f[f'particles/{particle_id}/trajectory'][:]
                trajectory_positions = trajectory[:, 1:dimensions+1].transpose()
                trajectory_velocities = trajectory[:, dimensions+1:dimensions*2+1].transpose()
                trajectories.append(trajectory_positions)
                if return_velocities:
                    velocities.append(trajectory_velocities)

        if return_velocities:
            return trajectories, velocities
        else:
            return trajectories

    def _plot_traj_colored(self, trajs: np.array, vels: np.array, ax: plt.Axes, dimensions: int, cmap: str = 'viridis')\
            -> LineCollection:
        """
        Plots the trajectories as lines made of colored segments representing the magnitude of
        the corresponding velocity, and also adds a color bar.

        :param trajs: A list of np.arrays, each np.array describing the particle positions with a shape of (d, n_i),
            where d is the spatial dimensions and n_i the number of points along the trajectory of the i-th particle
        :param vels: A list of velocities, matching the trajs parameter (1<->1 correspondence).
        :param ax: The axis (some plt.Axes instance) to plot on.
        :param dimensions: The number of spatial dimensions.
        :param cmap: The colormap to use, 'viridis' by default.
        :return: The LineCollection that was plotted.
        """
        VMAG = None
        SEGMENTS = None

        for i in range(len(trajs)):
            t = trajs[i]
            v = vels[i]
            vmag = np.linalg.norm(v.transpose(), axis=1)

            if dimensions == 2:
                points = t.T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
            elif dimensions == 3:
                points = t.T.reshape(-1, 1, 3)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

            VMAG = vmag if VMAG is None else np.concatenate((VMAG, vmag), axis=0)
            SEGMENTS = segments if SEGMENTS is None else np.concatenate((SEGMENTS, segments), axis=0)

        if dimensions == 3:
            from mpl_toolkits.mplot3d.art3d import Line3DCollection
            lc = Line3DCollection(SEGMENTS, cmap=cmap)
        else:
            lc = LineCollection(SEGMENTS, cmap=cmap)

        lc.set_array(VMAG)
        ax.add_collection(lc)
        ax.autoscale_view(True, True)
        ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='both')

        cb = plt.colorbar(lc, ax=ax)
        cb.set_label("velocity magnitude in m/s")

        return lc

    def _plot_detector_hits(self, h5f: h5py.File, ax: plt.Axes, dimensions: int) -> None:
        """
        Plots the detector hit positions from an HDF5 results file.

        :param h5f: A h5py.File instance
        :param ax: The axis to plot on
        :param dimensions: The number of spatial dimensions
        :return: Nothing
        """
        if 'detector_hits' in h5f:
            for detector_id in h5f['detector_hits']:
                hits = h5f['detector_hits'][detector_id][:]
                if dimensions == 3:
                    ax.scatter(hits[:, 0], hits[:, 1], zs=hits[:, 2], s=10, color='yellow')
                elif dimensions == 2:
                    ax.scatter(hits[:, 0], hits[:, 1], s=10, color='yellow')

    def _plot_positions(self, h5f: h5py.File, ax: plt.Axes, dimensions: int) -> None:
        """
        Plots the initial and final particle positions from an HDF5 results file.

        :param h5f: A h5py.File instance
        :param ax: The axis to plot on
        :param dimensions: The number of spatial dimensions
        :return: Nothing
        """
        particles = [h5f['particles'][k] for k in h5f['particles']]

        initial_positions = np.array([p['initial_position'][:dimensions] for p in particles])
        final_positions = np.array([p['final_position'][:dimensions] for p in particles])
        initial_positions_kwargs = {}
        final_positions_kwargs = {}

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.scatter(initial_positions[:, 0], initial_positions[:, 1], s=3, c='black', **initial_positions_kwargs)
        # FIXME these are often buggy so let's not plot them for now
        ax.scatter(final_positions[:, 0], final_positions[:, 1], s=3, c='green', **final_positions_kwargs)

    def visualize(self) -> Tuple[Figure, np.ndarray]:
        """
        Visualizes the results as described in the class.
        :return: A 2-tuple of (a plt.Figure instance, the array of plt.Axes instances created to plot)
        """
        with h5py.File(self.filename, 'r') as h5f:
            dimensions = int(h5f.attrs['dimensions'])
            if dimensions not in [2, 3]:
                raise ValueError("Can only plot 2D/3D!")

            # Construct the figure and axis
            fig = plt.figure()
            if dimensions == 3:
                # this import is not 'unused', it's a magic import that makes all the 3D stuff below possible
                from mpl_toolkits.mplot3d import Axes3D
                ax = fig.add_subplot(111, projection='3d')
                ax.set_zlabel('z')
            else:
                ax = fig.add_subplot(111)

            # Plot all trajectories that are present
            trajs, vels = self._get_trajectories(h5f, True, dimensions)
            self._plot_traj_colored(trajs, vels, ax, dimensions, cmap='viridis')

            # Plot all particle initial and final positions
            self._plot_positions(h5f, ax, dimensions)

            # Plot detector hits
            self._plot_detector_hits(h5f, ax, dimensions)

            T = np.concatenate(trajs, axis=1)
            X = T[0]
            Y = T[1]

            ax.set_xlim((np.min(X), np.max(X)))
            ax.set_ylim((np.min(Y), np.max(Y)))
            if dimensions == 3:
                Z = T[2]
                ax.set_zlim((np.min(Z), np.max(Z)))

        fig.tight_layout()
        return fig, np.array([ax]).reshape((1, 1))  # reshape to adhere to 2D axes ndarray


"""
### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
"""