import numpy as np
import scipy.optimize as opt
import emcee
from tqdm import tqdm

def log_likelihood(pred, obs, obs_err):
    return -0.5 * np.nansum((obs - pred)**2 / obs_err**2 + np.log(obs_err**2))

def neg_log_likelihood(pred, obs, obs_err):
    return 0.5 * np.nansum((obs - pred)**2 / obs_err**2 + np.log(obs_err**2))


# functions for running mcmc on a simple polynomial
def poly_log_likelihood(p, x, obs, sigma=1):
    """"
    Log-likeihood function for an Nth order polynomial.

    Parameters
    ----------
    p : array-like
        Containing the parameters to the polynomial, in the form [PN, P{N-1}, ..., P0].
        Passed directly to np.polyval().
    x : array-like
        The independent variable.
    obs : array-like
        Observed (dependent) data.
    sigma : array-like
        Standard deviations of the observed data.

    Returns
    -------
    float : log likelihood of function
    """
    pred = np.polyval(p, x)
    return log_likelihood(pred, obs, sigma)

def mcmc_poly(x, obs, sigma=1, order=1, nwalkers=32, niter=5000, start_sd=1e-2):
    """
    Run an MCMC sampler for a simple polynomial.

    Parameters
    ----------
    x : array-like
        The independent variable.
    obs : array-like
        Observed (dependent) data.
    sigma : array-like
        Standard deviations of the observed data.
    order : int
        The order of the polynomial to evaluate.
    nwalkers : int
        The number of walkers
    niter : int
        The number of iterations
    start_sd : float
        The spread in the initial conditions. Function minimum us multiplied
        by normally distributed random numbers with a mean of 1, and a
        standard deviation of start_sd.

    Returns
    -------
    emcee.ensemble.EnsembleSampler : MCMC sampler object
    """
    p0 = [0] * order
    
    negloglik = lambda p, x, obs, sigma: -poly_log_likelihood(p, x, obs, sigma)
    
    fmin = opt.minimize(negloglik, p0, args=(x, obs, sigma))
    start = np.random.normal(1, start_sd, (nwalkers, 1)) * fmin.x
    
    sampler = emcee.EnsembleSampler(nwalkers, order, poly_log_likelihood, args=(x, obs, sigma))
    
    for i in tqdm(sampler.sample(p0=start, iterations=niter), total=niter, desc='Sampling'):
        pass
    
    return sampler

def mcmc_fn(x, obs, loglik, p0, sigma=1, nwalkers=32, niter=5000, start_sd=1e-2):
    """
    Run an MCMC sampler for a simple polynomial.

    Parameters
    ----------
    x : array-like
        The independent variable.
    obs : array-like
        Observed (dependent) data.
    sigma : array-like
        Standard deviations of the observed data.
    loglik : int
        A function for form f(p, x, obs, sigma), which returns the log likelihood.
    p0 : array-like
        Initial geuss paramters for function minimisation
    nwalkers : int
        The number of walkers
    niter : int
        The number of iterations
    start_sd : float
        The spread in the initial conditions. Function minimum us multiplied
        by normally distributed random numbers with a mean of 1, and a
        standard deviation of start_sd.

    Returns
    -------
    emcee.ensemble.EnsembleSampler : MCMC sampler object
    """
    ndim = len(p0)
    
    negloglik = lambda p, x, obs, sigma: -loglik(p, x, obs, sigma)
    
    fmin = opt.minimize(negloglik, p0, args=(x, obs, sigma))
    start = np.random.normal(1, start_sd, (nwalkers, 1)) * fmin.x
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, loglik, args=(x, obs, sigma))
    
    for i in tqdm(sampler.sample(p0=start, iterations=niter), total=niter, desc='Sampling'):
        pass
    
    return sampler

def mcmc_pred(fn, x, p):
    """
    Predict values of function from MCMC chains.

    Parameters
    ----------
    fn : function
        With form f(p, x)
    x : np.ndarray
        A 1D array containing x values to evaluate the function over.
    p : array-like
        An array of shape (N, len(p)), where N is the number of MCMC
        samples.

    Returns
    -------
    np.ndarray : Containing evaluated function, of shape (len(x), N)
    """
    return fn(p.T, x[:, np.newaxis])

def mcmc_pval_pred(x, fchain):
    return np.polyval(fchain.T, x[:, np.newaxis])

def mcmc_CI(a, CI=0.95, axis=1):
    """
    Calculate median and confidence intervals from MCMC prediction.

    Parameters
    ----------
    a : np.ndarray
        Array containing output of mcmc_pred, of shape (len(x), N),
        where N is the number of MCMC samples.
    CI : float
        The confidence interval to report.
    
    Returns
    -------
    tuple : containing arrays of (lower, median, upper) of shape (len(x),)
    """
    if CI > 1:
        CI /= 100.
    return np.quantile(a, [(1 - CI) / 2, .50, 1 - (1 - CI) / 2], axis=axis)