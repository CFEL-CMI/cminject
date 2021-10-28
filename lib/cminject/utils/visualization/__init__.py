#!/usr/bin/env python3
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMInject
#
# Copyright (C) 2018,2020 CFEL Controlled Molecule Imaging group
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# If you use this program for scientific work, you should correctly reference it; see the LICENSE.md file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Code for visualization of virtual experiment results.
"""

import math
import warnings
from itertools import repeat
from typing import Optional, List, Union, Sequence, Dict, Callable

import numpy as np

from matplotlib import pyplot as plt, patheffects
from matplotlib.animation import Animation, FuncAnimation
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D

from cminject.utils.args import parse_multiple_dimensions_description, DimensionDescription


def _get_axis(ax: Optional[plt.Axes], dims: int):
    """
    Returns an appropriate plt.Axes instance for the given number of dimensions (``dims``).
    :param ax: (Optional) An existing axis. If it is incompatible with the number of dimensions,
      a warning will be printed. In any case, this axis will be returned unchanged.
    :param dims: The number of dimensions.
    :return: An appropriate plt.Axes instance for the given number of dimensions, if possible.
    """
    if ax is None:
        if dims not in [2, 3]:
            warnings.warn(f"I don't know what kind of axis to create for a {dims}D plot... Will create a 2D axis.")

        fig = plt.figure()
        if dims == 3:
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = fig.gca()
    else:
        if dims == 3 and not isinstance(ax, Axes3D):
            warnings.warn("You seem to have passed a 2D axis for 3D data! The plot may look incorrect.")
    return ax


def plot_trajectories(trajectories: Sequence[np.array], ax: Optional[plt.Axes] = None, **plot_kwargs):
    """
    Plots multiple trajectories as simple line plots (plt.plot).

    :param trajectories: An iterable of trajectories.
    :param ax: The axis to plot on; if None, a new figure and axis will be created.
    :param plot_kwargs: Keyword arguments that will be passed directly to ax.plot.
    :return: A list of the values returned by each ax.plot call.
    """
    if trajectories is None or not any(True for _ in trajectories):
        warnings.warn('No trajectories to plot were passed!')
        return None

    sample_traj = trajectories[0]

    julia_mode = 'position' not in sample_traj.dtype.fields
    if julia_mode:
        dims = 2
        ax = _get_axis(ax, dims)
    else:  # assume 2-D x/z from Julia code for now...
        dims = sample_traj['position'].shape[1]
        ax = _get_axis(ax, dims)

    plots = []
    for traj in trajectories:
        show = traj['position'].T if not julia_mode else (traj['x'], traj['z'])
        plots.append(ax.plot(*show, **plot_kwargs))
    return plots


def plot_trajectories_colored(trajectories: Sequence[np.array], ax: Optional[plt.Axes] = None,
                              autoscale: bool = True, **line_collection_kwargs):
    """
    Plots multiple trajectories, segments colored by the magnitude of the (approximate) velocities.

    :param trajectories: A sequence of trajectories.
    :param ax: The axis to plot on; if None, a new figure and axis will be created.
    :param autoscale: If True, the axes will be scaled to fit the created LineCollection/Line3DCollection.
    :param line_collection_kwargs: Keyword arguments that will be passed directly to the Line[3D]Collection constructor.
    :return: The Line[3D]Collection object representing all plotted lines.
    """
    if trajectories is None or not any(True for _ in trajectories):
        warnings.warn('No trajectories to plot were passed!')
        return None
    dims = trajectories[0]['position'].shape[1]
    assert dims == 2 or dims == 3, f"Don't know how to plot a {dims}-dimensional trajectory!"
    ax = _get_axis(ax, dims)

    all_vmags = None
    all_segments = None
    for traj in trajectories:
        pos = traj['position']
        t = traj['time']

        diff = np.diff(pos, axis=0)
        dt = np.diff(t)
        vmag = np.where(dt == 0, 0, np.divide(np.linalg.norm(diff, axis=1), dt))

        if dims == 2:
            points = pos.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
        elif dims == 3:
            points = pos.reshape(-1, 1, 3)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

        all_vmags = vmag if all_vmags is None else np.concatenate((all_vmags, vmag), axis=0)
        # noinspection PyUnboundLocalVariable
        all_segments = segments if all_segments is None else np.concatenate((all_segments, segments), axis=0)

    if dims == 3:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        lc = Line3DCollection(all_segments, **line_collection_kwargs)
    else:
        lc = LineCollection(all_segments, **line_collection_kwargs)

    lc.set_array(all_vmags)
    ax.add_collection(lc)

    if autoscale:
        # We could use ax.autoscale_view() in theory, but, alas: https://github.com/matplotlib/matplotlib/issues/17130
        mins = np.min(all_segments, axis=(0, 1))
        maxs = np.max(all_segments, axis=(0, 1))
        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        if dims == 3:
            ax.set_zlim(mins[2], maxs[2])

    return lc


def plot_trajectories_movie(trajectories: Sequence[np.array], dimension_description: DimensionDescription = "x,z",
                            delay: int = 100, n_previous: int = 10, ax: plt.Axes = None,
                            animation_kwargs: Dict = None, plot_kwargs: Dict = None,
                            update_fn: Callable[[int, List[plt.Artist]], List[plt.Artist]] = None)\
        -> (Animation, plt.Figure, int):
    """
    Creates a trajectory movie via a :class:`matplotlib.animation.Animation`. Plots each trajectory as a curve
    of the last ``n_previous`` positions.

    :param trajectories: A sequence of trajectories.
    :param dimension_description: A dimension descriptionf of a pair or triple of dimensions.
    :param delay: The delay between each animation frame.
    :param n_previous: The number of previous positions to include in each curve. When equal to one, each curve
    :param ax: (optional) a specific Axes instance to plot on. If not passed, a new Figure and Axes are created.
    :param animation_kwargs: (optional) arguments to pass to :class:`matplotlib.animation.FuncAnimation`
    :param plot_kwargs: (optional) arguments to pass to the initial plt.plot() calls that draw each line
    :param update_fn: (optional) An update function to be called as the last step in each animation iteration.
      Receives the current iteration number and the list of artists that the default update function returns.
      Must return a list of artists, just like animation update functions.
    :return: A 3-tuple of (created Animation instance, created Figure instance, number of animation frames).
    """
    extractor, n = parse_multiple_dimensions_description(dimension_description, n=(2, 3))\
        if type(dimension_description) is str else dimension_description
    positions = [extractor(traj) for traj in trajectories]
    maxlen = max(len(traj) for traj in trajectories)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111) if n == 2 else fig.add_subplot(111, projection='3d')
    else:
        fig = ax.get_figure()

    empty_plot_arg = [[]] * n
    lines = [plt.plot(*empty_plot_arg, **(plot_kwargs or {}))[0] for _ in range(len(trajectories))]
    labels = dimension_description.split(',')
    for label_setter, label in zip(('set_xlabel', 'set_ylabel', 'set_zlabel'), labels):
        getattr(ax, label_setter)(label)

    update_fn = update_fn or (lambda _, artists: artists)  # if update_fn is not supplied: just return artists

    def _init():
        for lim_setter, idx in zip(('set_xlim', 'set_ylim', 'set_zlim'), range(n)):
            getattr(ax, lim_setter)(min(np.min(pos[idx]) for pos in positions),
                                    max(np.max(pos[idx]) for pos in positions))
        return lines

    def _update(i):
        for k, ln in enumerate(lines):
            lower = max(0, i - n_previous - 1)
            ln.set_data(positions[k][0][lower:i], positions[k][1][lower:i])
            if n == 3:
                ln.set_3d_properties(positions[k][2][lower:i])
        return update_fn(i, lines)

    animation_kwargs = animation_kwargs or {}
    animation_kwargs = {'frames': maxlen, 'interval': delay, 'repeat': True, **animation_kwargs}
    ani = FuncAnimation(fig, _update, init_func=_init, **animation_kwargs)
    return ani, fig, maxlen


def plot_detector(detector: np.array, dimension_description: DimensionDescription, ax: plt.Axes, **kwargs):
    """
    Plots 1D/2D histograms of one or a pair of measured quantities at a single detector.

    :param detector: An np.ndarray representing all measurements (hits) made at the detector.
    :param dimension_description: A string like :func:`cminject.utils.args.parse_multiple_dimensions_description`
      accepts or a Callable as :func:`cminject.utils.args.parse_multiple_dimensions_description` returns.
      If it represents a pair of dimensions, a 2D histogram will be plotted, otherwise a 1D histogram will be plotted.
    :param ax: The axis to plot on.
    :param kwargs: Keyword args that will be passed directly to the plt.hist/hist2d call.
    :return: The result of the plt.hist/hist2d call made to plot the 1D/2D histogram.
    """
    extractor, n = parse_multiple_dimensions_description(dimension_description, n=(1, 2))\
        if type(dimension_description) is str else dimension_description
    result = extractor(detector)

    auto_bins = min(max(int(len(detector) / 100), 5), 250)
    if n == 2:
        kwargs2d = {'bins': auto_bins, 'cmin': 1, 'norm': LogNorm(), **kwargs}
        return ax.hist2d(result[0], result[1], **kwargs2d)
    else:  # n == 1
        kwargs1d = {'bins': auto_bins, **kwargs}
        return ax.hist(result, **kwargs1d)


def plot_detectors(detectors: Dict[str, np.array], dimension_description: str,
                   axes: Union[None, plt.Axes, List[plt.Axes]] = None, **kwargs):
    """
    Plots 1D/2D histograms of one or a pair of measured quantities at multiple detectors.

    :param detectors: The list of detectors, i.e., a list of what :func:`plot_detector` accepts.
    :param dimension_description: See docs for :func:`plot_detector`.
    :param axes: A list of axes or a single axis (plt.Axes instance). If one axis is passed, all plots will be made
      on that same axis; if a list of axes is passed, each detector will be plotted on the corresponding axis (by list
      position). This also means that if multiple axes are passed, len(axes) must be equal to len(detectors).
    :param kwargs: Keyword args that will be passed directly to the plt.hist/hist2d calls.
    :return: A list of all return values of :func:`plot_detector`.
    """
    n = len(detectors)
    if n == 0:
        return []

    if axes is None:
        xdim, ydim = round(math.sqrt(n)), math.ceil(math.sqrt(n))
        fig, axes = plt.subplots(xdim, ydim, sharex=True, sharey=True)
        axes = axes.flatten() if isinstance(axes, np.ndarray) else repeat(axes)
        fig.suptitle(dimension_description)
        if not isinstance(axes, repeat) and len(axes) > n:
            for ax in axes[n:]:
                ax.remove()
    else:
        assert len(axes) == len(detectors) or type(axes) is plt.Axes, \
            "Number of axes must match number of detectors or there must only be one axis given"
        if type(axes) is plt.Axes:
            axes = repeat(axes)

    extractor = parse_multiple_dimensions_description(dimension_description, n=(1, 2))
    plots = [plot_detector(detectors[d_id], extractor, ax, **kwargs) for ax, d_id in zip(axes, detectors.keys())]
    for ax, d_id in zip(axes, detectors.keys()):
        ax.axhline(0, color='k', alpha=.3, linewidth=1)
        ax.axvline(0, color='k', alpha=.3, linewidth=1)

        txt = ax.text(0.5, 0.87, d_id, transform=ax.transAxes, fontdict={'size': 9, 'color': 'black'}, ha='center')
        txt.set_path_effects([patheffects.withStroke(linewidth=3, foreground=(1,1,1,0.7))])
        ax.ticklabel_format(scilimits=(-2, 2))

    return plots
