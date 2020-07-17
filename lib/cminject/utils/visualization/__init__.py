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

import warnings
from itertools import repeat
from typing import Optional, Iterable, Callable, Any, List, Union
import math

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D

from cminject.utils.args import parse_dimension_description


def _get_axis(ax: Optional[plt.Axes], dims: int):
    if ax is None:
        fig = plt.figure()
        if dims == 3:
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = fig.gca()
    else:
        if dims == 3 and not isinstance(ax, Axes3D):
            warnings.warn("You seem to have passed a 2D axis for 3D data! The plot may look incorrect.")
    return ax


def plot_trajectories(trajectories: Iterable[np.array], ax: Optional[plt.Axes] = None, **plot_kwargs):
    if not any(True for _ in trajectories):
        return None
    dims = trajectories[0]['position'].shape[1]
    ax = _get_axis(ax, dims)

    plots = []
    for traj in trajectories:
        pos = traj['position']
        plots.append(ax.plot(*pos.T, **plot_kwargs))
    return plots


def plot_trajectories_colored(trajectories: Iterable[np.array], ax: Optional[plt.Axes] = None,
                              **line_collection_kwargs):
    if not any(True for _ in trajectories):
        return None
    dims = trajectories[0]['position'].shape[1]
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
        all_segments = segments if all_segments is None else np.concatenate((all_segments, segments), axis=0)

    if dims == 3:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        lc = Line3DCollection(all_segments, **line_collection_kwargs)
    else:
        lc = LineCollection(all_segments, **line_collection_kwargs)

    lc.set_array(all_vmags)
    ax.add_collection(lc)
    ax.autoscale_view(True, True)
    return lc


def plot_detector(detector: np.array, dimension_description: Union[str, Callable[[np.array], Any]],
                  ax: plt.Axes, **kwargs):
    if type(dimension_description) is str:
        extractor = parse_dimension_description(dimension_description)
    else:
        extractor = dimension_description

    result = extractor(detector)
    if type(result) is tuple:
        ax.hist2d(result[0], result[1], **kwargs)
    else:
        ax.hist(result, **kwargs)


def plot_detectors(detectors: List[np.array], dimension_description: str,
                   axes: Union[None, plt.Axes, List[plt.Axes]] = None, **kwargs):
    if axes is None:
        n = len(detectors)
        xdim, ydim = round(math.sqrt(n)), math.ceil(math.sqrt(n))
        fig, axes = plt.subplots(xdim, ydim)
        axes = axes.flatten() if isinstance(axes, np.ndarray) else repeat(axes)
        fig.suptitle(dimension_description)
        if len(axes) > n:
            for ax in axes[n:]:
                ax.remove()
    else:
        assert len(axes) == len(detectors) or type(axes) is plt.Axes, \
            "Number of axes must match number of detectors or there must only be one axis given"
        if type(axes) is plt.Axes:
            axes = repeat(axes)

    extractor = parse_dimension_description(dimension_description)
    for ax, detector in zip(axes, detectors):
        plot_detector(detector, extractor, ax, **kwargs)


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
