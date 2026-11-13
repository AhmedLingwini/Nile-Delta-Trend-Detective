"""Trend tools used by combine_power.py.

mann_kendall   Theil-Sen slope + Mann-Kendall test, p-value adjusted for lag-1 autocorrelation
bh_qvalues     Benjamini-Hochberg FDR q-values
season_series  seasonal means from a (points, years, 12) array of monthly values
"""
import numpy as np
from scipy.stats import norm, theilslopes


def mann_kendall(t, y, maxlag=8):
    """Return dict(slope, p_raw, p, r1, n) or None when fewer than 8 valid points.

    slope  Theil-Sen slope, in units of y per unit of t
    p_raw  Mann-Kendall two-sided p-value (no autocorrelation correction)
    p      p-value after inflating the variance of S by (1 + r1) / (1 - r1)
           when the lag-1 autocorrelation r1 of the detrended series is positive
    r1     lag-1 autocorrelation of the series after removing the Theil-Sen trend
    maxlag kept for compatibility with older calls; only lag 1 is used
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    n = len(y)
    if n < 8:
        return None
    slope = theilslopes(y, t)[0]
    s = np.sign(y[None, :] - y[:, None])[np.triu_indices(n, 1)].sum()
    var = n * (n - 1) * (2 * n + 5) / 18.0

    def pval(v):
        z = (s - np.sign(s)) / np.sqrt(v) if s != 0 else 0.0
        return 2 * (1 - norm.cdf(abs(z)))

    d = y - slope * t
    d = d - d.mean()
    r1 = float((d[:-1] * d[1:]).sum() / (d ** 2).sum())
    infl = max(1.0, (1 + r1) / (1 - r1))
    return dict(slope=float(slope), p_raw=float(pval(var)), p=float(pval(var * infl)), r1=r1, n=n)


def bh_qvalues(p):
    """Benjamini-Hochberg q-values for an array of p-values."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank, i in zip(range(m, 0, -1), order[::-1]):
        prev = min(prev, p[i] * m / rank)
        q[i] = prev
    return q


def season_series(arr, months):
    """arr: (points, years, 12) monthly values -> (points, years) mean over `months`."""
    with np.errstate(all="ignore"):
        return np.nanmean(np.asarray(arr, float)[:, :, months], axis=2)
