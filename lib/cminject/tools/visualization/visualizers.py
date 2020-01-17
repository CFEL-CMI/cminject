#!/usr/bin/env python3
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
from typing import List, Tuple, Any, Callable, Union, Dict

import h5py
import numpy as np
from cminject.definitions.result_storages.hdf5 import HDF5ResultStorage
from matplotlib import patheffects
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection
from matplotlib.colorbar import Colorbar
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.art3d import Line3D


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
    def visualize(self, *args, **kwargs) -> Any:
        """
        A highly generic method to construct a visualization based on the filename passed at construction.

        Note that it's not conceptually required that calling this method shows anything on screen, see e.g.
        the docstring for MatplotlibVisualizer as to why this is.

        :return: Any set of objects ne. See e.g. MatplotlibVisualizer for a more concretely types
        """
        pass


DimensionType = Union[str, Callable[[List[np.array]], np.array]]


def _get_dims(description: List[str], dims: List[DimensionType], hits: List[np.array]) -> List[np.array]:
    return [
        dim(hits) if callable(dim)
        else hits[np.where(description == dim)[0][0]]
        for dim in dims
    ]


class DetectorHistogramVisualizer(Visualizer):
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

    def __init__(self, filename: str, dimension_pairs: List[Tuple[DimensionType, DimensionType]],
                 colormap: str = 'viridis', highlight_origin: bool = True,
                 bins: int = 100, bins_1d: int = 30, use_tex: bool = False,
                 hist1d_kwargs: Dict[Any, Any] = None, hist2d_kwargs: Dict[Any, Any] = None):
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
        :param use_tex: Whether to use TeX for rendering text.
        """
        super().__init__(filename)
        self.dimension_pairs = dimension_pairs
        self.colormap = colormap
        self.highlight_origin = highlight_origin
        self.bins = bins
        self.bins_1d = bins_1d
        self.use_tex = use_tex
        self.hist1d_kwargs = hist1d_kwargs or {}
        self.hist2d_kwargs = hist2d_kwargs or {}

    @staticmethod
    def vis1d(bins, ax, y, **hist_kwargs) -> None:
        h = ax.hist(y, bins=bins, **hist_kwargs)
        return h

    @staticmethod
    def vis2d(bins, fig, ax, x, y, **hist2d_kwargs):
        h2d = ax.hist2d(x, y, bins=bins, **hist2d_kwargs)
        img = h2d[3]
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cax.ticklabel_format(scilimits=(-4, 4), useMathText=True)
        fig.colorbar(img, cax=cax)
        return h2d

    def visualize(self) -> Tuple[Figure, np.ndarray]:
        """
        Visualize the results as a plot of (potentially multiple) 1D/2D histograms.
        There will be m*n subplots arranged on a grid, with m being the number of detectors in the file,
        and n being the number of dimension pairs given at construction of this instance.

        :return: A Matplotlib figure and a 2-dimensional ndarray of all axes present in the figure.
        """
        if self.use_tex:
            plt.rc('text', usetex=True)
            plt.rc('axes', unicode_minus=False)

        with h5py.File(self.filename, 'r') as h5f:
            if 'detector_hits' not in h5f:
                fig = plt.figure()
                return fig, fig.gca()

            detectors = h5f['detector_hits']
            detector_count = len(detectors)
            detector_keys = list(sorted(detectors.keys(), key=int))
            dimpair_count = len(self.dimension_pairs)

            fig, axes = plt.subplots(
                nrows=detector_count, ncols=dimpair_count,
                figsize=(dimpair_count, detector_count),
                gridspec_kw={'hspace': .7, 'wspace': .7}
            )
            axes = np.array(axes).reshape((detector_count, dimpair_count))
            # Plot a hist2d in each subplot in the grid :)
            for j in range(dimpair_count):
                for i in range(detector_count):
                    # Get all hits for this detector
                    hits = detectors[detector_keys[i]]
                    description = hits.attrs['description']
                    hits = hits[:].transpose()

                    # Find the dimensions' indices and get the correct axis to plot on
                    dim_a, dim_b = self.dimension_pairs[j]
                    try:
                        x, y = _get_dims(description, [dim_a, dim_b], hits=hits)
                    except IndexError as e:
                        raise ValueError(f"Could either not find dimension {dim_a} or dimension {dim_b}!")

                    ax = axes[i, j]
                    ax.ticklabel_format(scilimits=(-2, 2))
                    ax.ticklabel_format(scilimits=(-2, 2))
                    ax.autoscale(enable=True, tight=False)  # enable auto-scaling here :)

                    # Set the title of the column and decide on the type of histogram
                    if i == 0:
                        title = f"{dim_a} / {dim_b}"
                        title = f"${title}$" if self.use_tex else title
                        ax.set_title(title)
                    if j == 0:
                        ax.set_ylabel(f"Det. {detector_keys[i]}", rotation=90, size='large')

                    one_dimensional = False
                    if np.all(x == x[0]):  # X values are all equal, plot 1D histogram of Y
                        one_dimensional = True
                        self.vis1d(self.bins_1d, ax, y, **self.hist1d_kwargs)
                    elif np.all(y == y[0]):  # Y values are all equal, plot 1D histogram of X
                        one_dimensional = True
                        self.vis1d(self.bins_1d, ax, x, **self.hist1d_kwargs)
                    else:
                        self.vis2d(self.bins, fig, ax, x, y)

                    # Origin highlighting for 2D histograms
                    if self.highlight_origin and not one_dimensional:
                        # Disable auto-scaling here: If 0.0,0.0 is outside the existing bounds, that's fine, we don't
                        # want to stretch the existing plot just to show the origin marker
                        ax.autoscale(enable=False)
                        ax.scatter([0.0], [0.0], s=10, c=[(1.0, 0.0, 0.0, 0.5)], marker='x')
        fig.tight_layout()
        return fig, axes


class TrajectoryVisualizer(Visualizer):
    """
    Visualizes the initial and final positions of all particles stored in an HDF5 file, along with their trajectories,
    in a 2D or 3D plot (based on the dimensionality of the stored particle positions). Detector hits are also shown.

    Useful to get a general idea of where the particles started and went, but not very useful for data analysis:
    the plots can get rather chaotic and slow (especially in 3D).
    """
    def __init__(self, *args, n_samples: int = None,
                 fig: plt.Figure = None, ax: plt.Axes = None,
                 colored: bool = False, colorbar: bool = True,
                 make_x_absolute: bool = False,
                 traj_kwargs: Dict[Any, Any] = None,
                 scatter: bool = False, scatter_kwargs: Dict[Any, Any] = None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.n_samples = n_samples
        self.fig = fig
        self.ax = ax
        self.colorbar = colorbar

        self.colored = colored
        self.make_x_absolute = make_x_absolute
        self.traj_kwargs = traj_kwargs or {}
        self.scatter = scatter
        self.scatter_kwargs = scatter_kwargs or {}

    def visualize(self) -> Tuple[Figure, Axes]:
        """
        Visualizes the results as described in the class.
        :return: A 2-tuple of (the plt.Figure instance, the plt.Axes instance created to plot)
        """
        storage = HDF5ResultStorage(self.filename, mode='r')
        dimensions = storage.get_dimensions()

        if dimensions not in [2, 3]:
            raise ValueError("Can only plot 2D/3D!")

        # Construct the figure and axis
        ax = self.ax
        fig = self.fig
        if fig is None:
            fig = plt.figure()
        if ax is None:
            if dimensions == 3:
                # this import is not 'unused', it's a magic import that makes all the 3D stuff below possible
                # noinspection PyUnresolvedReferences
                from mpl_toolkits.mplot3d import Axes3D
                ax = fig.add_subplot(111, projection='3d')
                ax.set_zlabel('z')
            else:
                ax = fig.add_subplot(111)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        # Plot all trajectories that are present
        print("Getting trajectories...")
        trajectories = storage.get_trajectories_iterator(self.n_samples)

        print("Showing trajectories...")
        if trajectories:
            if self.colored:
                self.plot_traj_colored(
                    trajectories, ax, dimensions, colorbar=self.colorbar, **self.traj_kwargs
                )
            else:
                self.plot_traj_monochrome(
                    trajectories, ax, dimensions, **self.traj_kwargs
                )

        # Plot all particle initial and final positions
        if self.n_samples is None and self.scatter:
            print("Showing initial and final positions...")
            initial, final = storage.get_initial_and_final_positions()
            self.plot_positions(initial, final, ax, dimensions, **self.scatter_kwargs)

        # Plot detector hits
        if self.scatter:
            print("Showing detected positions...")
            detectors = storage.get_detectors()
            self.plot_detector_hits(detectors, ax, dimensions, **self.scatter_kwargs)

        ax.autoscale_view(True, True)
        fig.tight_layout()
        return fig, ax

    @staticmethod
    def plot_traj_monochrome(trajectories: np.array, ax: plt.Axes, dimensions: int, **plot_kwargs)\
            -> List[Union[Line2D, Line3D]]:
        """
        Plots the trajectories as monochromatic lines (curves).

        :param trajectories: A list of np.arrays, each np.array containing the full trajectory of the corresponding
            particle.
        :param ax: The axis (some plt.Axes instance) to plot on.
        :param dimensions: The number of spatial dimensions.
        :param plot_kwargs: The kwargs to pass to plt.plot.
        :return: A list of line instances (Line2D or Line3D) that were plotted.
        """
        lines = []

        for traj in trajectories:
            lines.append(ax.plot(*traj[1:dimensions+1], **plot_kwargs))

        ax.autoscale_view(True, True)
        ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='both')

        return lines

    @staticmethod
    def plot_traj_colored(trajectories: np.array, ax: plt.Axes, dimensions: int, colorbar: bool = True, **line_kwargs)\
            -> Tuple[LineCollection, Colorbar]:
        """
        Plots the trajectories as lines made of colored segments representing the magnitude of
        the corresponding velocity, and also adds a color bar.

        :param trajectories: A list of np.arrays, each np.array containing the full trajectory of the corresponding
            particle.
        :param ax: The axis (some plt.Axes instance) to plot on.
        :param dimensions: The number of spatial dimensions.
        :param colorbar: Whether to plot a color bar. True by default.
        :param line_kwargs: The kwargs to pass to LineCollection / Line3DCollection.
        :return: A 2-tuple of the LineCollection and the colorbar that was plotted. The returned colorbar may be None,
            if the colorbar argument passed to this function was False.
        """
        VMAG = None
        SEGMENTS = None

        for traj in trajectories:
            t = traj[0]
            p = traj[1:dimensions+1]

            diff = np.diff(p)
            dt = np.diff(t)
            vmag = np.where(dt == 0, 0, np.divide(np.linalg.norm(diff, axis=0), dt))

            if dimensions == 2:
                points = p.T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
            elif dimensions == 3:
                points = p.T.reshape(-1, 1, 3)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

            VMAG = vmag if VMAG is None else np.concatenate((VMAG, vmag), axis=0)
            SEGMENTS = segments if SEGMENTS is None else np.concatenate((SEGMENTS, segments), axis=0)

        if dimensions == 3:
            from mpl_toolkits.mplot3d.art3d import Line3DCollection
            lc = Line3DCollection(SEGMENTS, **line_kwargs)
        else:
            lc = LineCollection(SEGMENTS, **line_kwargs)

        lc.set_array(VMAG)
        ax.add_collection(lc)
        ax.autoscale_view(True, True)
        ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='both')

        if colorbar:
            cb = plt.colorbar(lc, ax=ax)
            cb.set_label("velocity magnitude in m/s")
        else:
            cb = None

        return lc, cb

    @staticmethod
    def plot_detector_hits(detectors: List[Tuple[Any, np.array]], ax: plt.Axes, dimensions: int, **scatter_kwargs)\
            -> None:
        """
        Plots the detector hit positions from an HDF5 results file.

        :param detectors: A list of detectors, each detector represented as a tuple of an identifier and an np.array
            of detected positions (+ other particle properties at time of detection)
        :param ax: The axis to plot on
        :param dimensions: The number of spatial dimensions
        :return: Nothing
        """
        if 'zorder' not in scatter_kwargs:
            scatter_kwargs.update({'zorder': 3})
        if 's' not in scatter_kwargs:
            scatter_kwargs.update({'s': 10})

        for i, hits in detectors:
            if dimensions == 3:
                ax.scatter(hits[0], hits[1], zs=hits[2], **scatter_kwargs)
            elif dimensions == 2:
                ax.scatter(hits[0], hits[1], **scatter_kwargs)

    @staticmethod
    def plot_positions(initial: np.array, final: np.array, ax: plt.Axes, dimensions: int, **scatter_kwargs) -> None:
        """
        Plots the initial and final particle positions from an HDF5 results file. The initial positions
        are plotted as a scatter plot of black points, the final positions as a scatter plot of green points.

        :param initial: An np.array containing all initial positions
        :param final: An np.array containing all final positions
        :param ax: The axis to plot on
        :param dimensions: The number of spatial dimensions
        :param scatter_kwargs: The kwargs to pass to plt.scatter.
        :return: Nothing
        """
        initial_kwargs, final_kwargs = scatter_kwargs.copy(), scatter_kwargs.copy()
        if dimensions == 3:
            initial_kwargs['zs'] = initial[2]
            final_kwargs['zs'] = final[2]

        ax.scatter(initial[0], initial[1], s=3, c='black', **initial_kwargs)
        ax.scatter(final[0], final[1], s=3, c='green', **final_kwargs)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
