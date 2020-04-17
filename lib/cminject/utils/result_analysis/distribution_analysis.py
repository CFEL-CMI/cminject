import warnings

import numpy as np
from scipy import signal, stats
from sklearn import linear_model


def get_x_percent_position(data, x):
    if len(data) == 0:
        return np.nan

    data = np.abs(data)
    counts, binedges = np.histogram(data, bins=np.linspace(0, np.max(data), 5000))
    total_counts = np.sum(counts)
    cumsum = np.cumsum(counts)
    # use the left bin edge as an acceptable approximation
    return binedges[np.argmax(cumsum > (x/100)*total_counts)]


def get_maxpeak_position(xs):
    # filter out outliers, which are most likely experimental noise
    # TODO find a decent cutoff for the Z-score
    zscores = np.abs(stats.zscore(xs))
    idxs = np.argwhere(zscores < 10).ravel()  # save the original indices for reconstruction
    xs = xs[idxs]  # use the indices for filtering

    # find peaks in the filtered data
    # TODO find good expected peak widths, or use find_peaks
    peaks = signal.find_peaks_cwt(xs, [10, 20, 30])
    if len(peaks) == 0:
        return np.nan

    # find the index of the maximum peak in the filtered data, then get the original index via idxs
    maxpeak_idx = np.argmax(xs[peaks])
    maxpeak = peaks[maxpeak_idx]
    maxpos = idxs[maxpeak]
    return maxpos


def _rotate(p, origin=(0, 0), angle=0.0):
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle), np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T - o.T) + o.T).T)


def autorotate_positions(positions):
    X, y = positions[:, 0].reshape(-1, 1), positions[:, 1]
    model = linear_model.LinearRegression()
    model.fit(X, y)
    angle = -np.arctan(model.coef_[0])
    positions_rot = _rotate(
        positions,
        origin=(0, model.intercept_),
        angle=angle
    )

    X_, y_ = positions_rot[:, 0].reshape(-1, 1), positions_rot[:, 1]
    model_rot = linear_model.LinearRegression()
    model_rot.fit(X_, y_)

    if model_rot.coef_[0] > 0.01:
        warnings.warn(f'Could not rotate positions well, returning original positions. '
                      f'Leftover ascent is {model_rot.coef_[0]}')
        return positions, (model, model_rot)

    return positions_rot, (model, model_rot)