"""
Misc plotting tools.
"""
import numpy as np

def rangecalc(x, pad=0.05):
    """
    Calculate padded axis limits based on x.abs

    Parameters
    ----------
    x : array-like
        The values to consider.
    pad : float
        The amount to pad the axis by, as a proportion of the data range.
    
    Returns
    -------
    (min, max) : tuple
    """
    mn = np.nanmin(x)
    mx = np.nanmax(x)
    rn = mx - mn
    return mn - pad * rn, mx + pad * rn


def spreadm(x, y, x_tol, y_tol, x_offset=None, y_offset=None, offset_mult=0.2):
    """
    spreadm reditributes a set of overlapping x/y points around their mean for display.

    Parameters
    ----------
    x, y : array-like
        The x and y arrays containing overlapping points
    x_tol, y_tol : float
        The overlap tolerance. Points are redistributed if any
        of them are closer to the mean than either tolerance.
    x_offset, y_offset : float
        Absolute x/y offsets to redistribute the points by.
        If None, offsets are calculated based on the number
        of points, and the sizes of the tolerances, following:
        offset = tolerance * offset_mult * n
        Defaults to None.
    offset_mult : float
        Used to automatically calculate displacement offsets. See above.

    Returns
    -------
    x_new, y_new, x_mean, y_mean

    If no points are outside the tolerances, or there is only one
    point, x_mean and y_mean are None.
    """
    x_mean = np.nanmean(x)
    y_mean = np.nanmean(y)

    x_diff = abs(x - x_mean)
    y_diff = abs(y - y_mean)

    if any(x_diff > x_tol) or any(y_diff > y_tol) or (len(x) >= 2):
        n = len(x)
        rad = 2 * np.pi / n

        rads = np.arange(n) * rad

        if x_offset is None:
            x_offset = x_tol * offset_mult * n
        if y_offset is None:
            y_offset = y_tol * offset_mult * n

        x_out = []
        y_out = []
        for r in rads:
            y_off = y_offset * np.cos(r)
            x_off = x_offset * np.sin(r)

            x_out.append(x_off)
            y_out.append(y_off)
        return x_mean + x_out, y_mean + y_out, x_mean, y_mean
    else:
        return x, y, None, None


def intervals(x, y, f, p, xn=None, interval_type='confidence', conflevel=0.95):
    """
    General function to calculate the confidence or
    prediction interval for a fitted dataset using
    equations laid out here:
    http://www.jerrydallal.com/LHSP/slr.htm

    Parameters:
        x, y: array-like
            raw data arrays
        f: function
            the fit function, of form f(x, *p)
        p: array-like
            the fit parameters
        xn: array-like
            the x-range over which to return the interval
            if None,
            returns a sequence over the original x range,
            with n=100.
        conflevel: float
            confidence level of prediction (default = 0.95)
        interval_type: string
            'confidence' returns confidence interval
            'prediction' returns prediction interval
    """
    import numpy as np
    from scipy.stats import t

    # set up calculated parameters
    alpha = 1. - conflevel  # significance level
    n = x.size  # data sample size

    if xn is None:
        xn = np.linspace(x.min(), x.max(), 100)

    # calculate predicted values
    yp = f(x, *p)  # on original x-scale
    ynp = f(xn, *p)  # on new x-scale

    # calculate standard error of estimate
    # (http://onlinestatbook.com/2/regression/accuracy.html)
    Se = (np.sum((y - yp)**2) / (n - 2))**.5

    # calculate quantile of Student's t distribution
    # for p = 1 - alpha/2 (UNSURE WHY!)
    q = t.ppf(1. - alpha / 2, n - 2)

    # distance from data centre, for formula see
    # see http://www.jerrydallal.com/LHSP/slr.htm
    dx = (xn - x.mean())**2 / np.sum((xn - xn.mean())**2)

    # calculate distance from prediction line
    if interval_type is 'confidence':
        dy = q * Se * np.sqrt(1 / n + dx)
    if interval_type is 'prediction':
        dy = q * Se * np.sqrt(1 + 1 / n + dx)

    return xn, ynp, ynp + dy, ynp - dy


def unitpicker(a, llim=0.1):
    """
    Calculate chemical units of a, such that no values are less than llim.
    """
    if isinstance(a, (np.ndarray, list)):
        a = np.nanmin(a)
    
    udict = {0: 'mol/mol',
             1: 'mmol/mol',
             2: '$\mu$mol/mol',
             3: 'nmol/mol',
             4: 'pmol/mol',
             5: 'fmol/mol'}
    a = abs(a)
    n = 0
    if a < llim:
        while a < llim:
            a *= 100
            n += 1
    return float(1000**n), udict[n]

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])