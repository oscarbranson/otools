"""
Peakshape library. Largely redundant, and present in other libraries.

e.g. scipy.stats, sklearn
"""
import numpy as np

def gaussian(x, area, cen, fwhm):
    """
    Gaussian Peak

    Parameters
    ----------
    x : array-like
    area : float
        Peak area
    cen : float
        Peak centre
    fwhm : float
        Peak full-width at half maximum

    Return
    ------
    y : array-like
    """
    return (area / fwhm * np.sqrt(4 * np.log(2) / np.pi) *
            np.exp(-4 * np.log(2) * ((x - cen) / fwhm)**2))

def gaussian_assym(x, area, cen, fwhm0, assym):
    """
    Gaussian Peak

    Parameters
    ----------
    x : array-like
    area : float
        Peak area
    cen : float
        Peak centre
    fwhm0 : float
        Peak full-width at half maximum for symmetric peak
    assym : float
        expoential coefficient defining variations of FWHM with distance from cen.
        0 = symmetric

    Return
    ------
    y : array-like
    """
    fwhm_x = (2 * fwhm0) / (1 + np.exp(assym * (x - cen)))

    return (area / fwhm_x * np.sqrt(4 * np.log(2) / np.pi) *
            np.exp(-4 * np.log(2) * ((x - cen) / fwhm_x)**2))


def lorentzian(x, area, cen, fwhm):
    """
    Lorentzian Peak

    Parameters
    ----------
    x : array-like
    area : float
        Peak area
    cen : float
        Peak centre
    fwhm : float
        Peak full-width at half maximum

    Return
    ------
    y : array-like
    """
    return (2 * area / (np.pi * fwhm)) / (1 + 4 * ((x - cen) / fwhm)**2)


def pvoigt(x, area, cen, fwhm, frac):
    """
    Pseudo-Voigt Peak

    Parameters
    ----------
    x : array-like
    area : float
        Peak area
    cen : float
        Peak centre
    fwhm : float
        Peak full-width at half maximum
    frac : float
        Proportion of Gaussian vs Lorentzian, where y = frac * L(x) + (1 - frac) * G(x)

    Return
    ------
    y : array-like
    """
    return frac * lorentzian(x, area, cen, fwhm) + (1 - frac) * gaussian(x, area, cen, fwhm)

# as a single function
def pvoigt_assym(x, area, cen, fwhm0, assym, frac):
    """
    Assymetric Pseudo-Voigt Peak

    Parameters
    ----------
    x : array-like
    area : float
        Peak area
    cen : float
        Peak centre
    fwhm0 : float
        Peak full-width at half maximum without skew
    assym : float
        expoential coefficient defining variations of FWHM with distance from cen.
    frac : float
        Proportion of Gaussian vs Lorentzian, where y = frac * L(x) + (1 - frac) * G(x)

    Return
    ------
    y : array-like
    """
    fwhm_x = (2 * fwhm0) / (1 + np.exp(assym * (x - cen)))
    return (frac * (2 * area / (np.pi * fwhm_x)) /
            (1 + 4 * ((x - cen) / fwhm_x)**2) +
            (1 - frac) * area / fwhm_x * np.sqrt(4 * np.log(2) / np.pi) *
            np.exp(-4 * np.log(2) * ((x - cen) / fwhm_x)**2))


# from https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/
def G(x, alpha):
    """ Return Gaussian line shape at x with HWHM alpha """
    return np.sqrt(np.log(2) / np.pi) / alpha\
                             * np.exp(-(x / alpha)**2 * np.log(2))

def L(x, gamma):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return gamma / np.pi / (x**2 + gamma**2)

def V(x, alpha, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)
