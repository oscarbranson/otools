"""
Tools for plotting (x,y) data and associated errors as 2D gaussians on an image.
"""
import numpy as np
from scipy import stats

def gauss2d(ps, xs, ys):
    """
    Generate a 2d gaussian on xs, ys

    Parameters
    ----------
    ps : (x, y, xe, ye) tuple
        Centre (x, y) and error (xe, ye) for points
    xs, ys : array-like
        Axes on which to plot data. Should be output of np.meshgrid.

    Returns
    -------
    d : array-like
        2d gaussian same size as xs and ys.
    """
    x, y, xe, ye = ps
    return stats.norm(x, xe).pdf(xs) * stats.norm(y, ye).pdf(ys)


def gaussplot(x, y, xe, ye, n=500, pad=0.1):
    """
    Plot seies of points with x and y errors as stacks of probability density functions.

    Parameters
    ----------
    x, y : array-like
        Data points to plot
    xe, ye : array-like
        Errors on data points
    n : int
        The size of the returned array (n, n)
    pad : float
        The proportion of the axis range to pad the axis.
    
    Returns
    -------
    xs, ys, d : tuple of array-like
        xs and ys are coordinate arrays, d is plot image.
    """
    xs, ys = np.meshgrid(np.linspace(x.min() - x.ptp() * pad,
                                     x.max() + x.ptp() * pad, n),
                         np.linspace(y.min() - y.ptp() * pad,
                                     y.max() + y.ptp() * pad, n))
    
    d = np.nansum(np.apply_along_axis(gauss2d, 0, np.vstack([x, y, xe, ye]), xs=xs, ys=ys), -1)
    
    return xs, ys, d