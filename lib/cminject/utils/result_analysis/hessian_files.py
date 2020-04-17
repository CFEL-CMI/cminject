import warnings

import h5py
import numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl

from .distribution_analysis import \
    get_x_percent_position, get_maxpeak_position, autorotate_positions


def _draw_boundary(ax, xspace, yspace):
    ax.axvline(np.min(xspace), linewidth=.5, color='black')
    ax.axvline(np.max(xspace), linewidth=.5, color='black')
    ax.axhline(np.min(yspace), linewidth=.5, color='black')
    ax.axhline(np.max(yspace), linewidth=.5, color='black')


class Slice:
    def __init__(self, points, left, right, xshift, zshift, zorigin=0.0):
        self.points = points
        self.left = left
        self.right = right
        self.xshift = xshift
        self.zshift = zshift
        self.zorigin = zorigin

    def get_x_percent_position(self, x):
        xdist = np.abs(self.get_xs(shift=True))
        return get_x_percent_position(xdist, x=x)

    def get_xs(self, shift=True):
        return self.points[1] - (self.xshift if shift else 0)

    def get_zs(self, shift=True):
        return self.points[0] - (self.zshift if shift else 0)

    def __str__(self):
        return f"<Slice@({self.left}, {self.right}), {self.points.shape[1]}>"

    def __repr__(self):
        return self.__str__()


class HessianFile:
    def __init__(self, filename, min_maxintensity=600,
                 resolution=(400, 400), extent=(740e-6, 740e-6), zorigin=0.0):
        self.filename = filename
        with h5py.File(filename, 'r') as f:
            self.integrated_intensities = f['IntegratedIntensity'][:].ravel()
            self.frame_numbers = f['FrameNumber'][:].ravel()
            self.maximum_intensities = f['MaximumIntensity'][:].ravel()
            # flip the data in z so we have the correct propagation direction
            self.positions = resolution[0] - f['Position'][:]  # assuming z is first dimension

        self.min_maxintensity = min_maxintensity
        self.extent = extent
        self.resolution_px = resolution
        self.resolution_m = tuple(extent[i] / resolution[i] for i in range(len(resolution)))
        self.xspace_px, self.yspace_px = [np.linspace(0, resolution[i], resolution[i] + 1) for i in
                                          range(2)]
        self.zorigin = zorigin

        psr, (model0, model1) = autorotate_positions(self.positions)
        psr_filter = psr[self.maximum_intensities > self.min_maxintensity]
        # Conceptual shift to x/z starting here
        hx, _ = np.histogram(psr_filter[:, 1], bins=self.yspace_px)
        hz, _ = np.histogram(psr_filter[:, 0], bins=self.xspace_px)
        self.hx = hx
        self.hz = hz
        self.xpeak, self.zpeak = [get_maxpeak_position(dat) for dat in (hx, hz)]
        self.rotated_positions = psr
        self.model0, self.model1 = model0, model1

    def px_to_m(self, px, axis=0, absolute=False):
        m = self.resolution_m[axis] * px
        if absolute and axis == 0:
            return m + self.zorigin
        else:
            return m

    def get_slices(self, n_slices):
        slices = []
        step = self.resolution_px[0] / n_slices  # along z

        for i in range(n_slices):
            left = i * step
            right = (i + 1) * step

            points = self.rotated_positions.T
            points = points[:, (left <= points[0]) & (points[0] < right)]
            slices.append(Slice(points, left, right,
                                xshift=self.xpeak, zshift=self.zpeak, zorigin=self.zorigin))

        return slices

    def plot_overview(self):
        ps = self.positions
        psr = self.rotated_positions

        # Plotting: Init
        fig, axs = plt.subplots(2, 2, figsize=(9, 9))
        axs = axs.ravel()

        # Plotting: Original -> Rotated
        axs[0].scatter(*ps.T, color='darkred', s=1, label='original')
        axs[0].scatter(*psr.T, color='C0', s=1, label='rotated')
        axs[0].plot(ps[:, 0], self.model0.predict(ps[:, 0].reshape((-1, 1))), color='black')
        axs[0].plot(psr[:, 0], self.model1.predict(psr[:, 0].reshape((-1, 1))), color='black')
        _draw_boundary(axs[0], self.xspace_px, self.yspace_px)

        axs[0].set_title('Original and auto-rotated scatterplot')
        axs[0].set_xlabel('z');
        axs[0].set_ylabel('x')
        axs[0].legend()

        # Plotting: Advanced overview
        scatter = axs[1].scatter(*psr.T, s=1, c=self.maximum_intensities,
                                 norm=mpl.colors.LogNorm(), alpha=.3)
        fig.colorbar(scatter, ax=axs[1])

        axs[1].axvline(self.zpeak, color='lime')
        axs[1].axhline(self.xpeak, color='lime')
        axs[1].bar(self.yspace_px[:-1], self.hz / np.max(self.hz) * np.max(self.yspace_px) / 4,
                   color='lime')
        axs[1].barh(self.xspace_px[:-1], self.hx / np.max(self.hx) * np.max(self.xspace_px) / 4,
                    color='lime')
        _draw_boundary(axs[1], self.xspace_px, self.yspace_px)
        axs[1].set_title('Intensities and intensity peaks along $z$, $x$')
        axs[1].set_xlabel('z');
        axs[1].set_ylabel('x')

        # Plotting: Peaks in x
        axs[2].hist(psr[:, 1], bins='auto', label='All', orientation='horizontal')
        axs[2].hist(psr[self.maximum_intensities > self.min_maxintensity][:, 1],
                    bins='auto', label='$I_{max}$ > 600', orientation='horizontal')
        axs[2].axhline(self.xpeak, color='lime')

        axs[2].set_title('x positions and peak')
        axs[2].set_ylabel('x');
        axs[2].set_xlabel('Occurrences')
        axs[2].legend()

        # Plotting: Peaks in z
        axs[3].hist(psr[:, 0], bins='auto', label='All')
        axs[3].hist(psr[self.maximum_intensities > self.min_maxintensity][:, 0],
                    bins='auto', label='$I_{max}$ > 600')
        axs[3].axvline(self.zpeak, color='lime')

        axs[3].set_title('z positions and peak')
        axs[3].set_xlabel('z');
        axs[3].set_ylabel('Occurrences')
        axs[3].legend()

        fig.tight_layout()

    def plot_focus_curve(self, n_slices=8, ax=None, percentage=70, zfac=1e3, xfac=1e6, **kwargs):
        slices = self.get_slices(n_slices=n_slices)
        slice_zs = np.array([(slice_.left + slice_.right) / 2 for slice_ in slices])
        slice_pNs = np.array([slice_.get_x_percent_position(x=percentage) for slice_ in slices])

        if ax is None:
            fig = plt.figure()
            ax = fig.gca()
        ax.set_xlabel(f'z distance [m/{zfac:.0g}]')
        ax.set_ylabel(f'x distance [m/{xfac:.0g}] from center until {percentage:3.0f}%')
        plot = ax.plot(
            (self.zorigin + (slice_zs - self.zpeak) * self.resolution_m[0]) * zfac,
            (slice_pNs * self.resolution_m[1]) * xfac,
            '-o',
            **kwargs
        )
        return plot