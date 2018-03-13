import numpy as np
from scipy import stats


def resample(N, *it):
    """
    Returns N samples from all values.

    Parameters
    ----------
    N : int or tuple of ints
        The shape of samples to return.
    *it : scipy.stats.rv_continuous objects or tuples
        If a distribution object, samples the distribution directly.
        If a tuple, takes the first value as the mean and the second
        as the standard deviation for sampling.
    
    Returns
    -------
    numpy.ndarray of shape (len(it), *N)
    """
    if isinstance(N, float):
        N = int(N)
    elif isinstance(N, (tuple, list, np.ndarray)):
        N = np.array(N, dtype=int)

    out = []
    for i in it:
        if isinstance(i, tuple):
            i = stats.norm(i[0], i[1])
        out.append(i.rvs(N))
    return np.array(out)