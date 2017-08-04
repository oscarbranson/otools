import numpy as np


def gaussian(x, area, cen, fwhm):
    return (area / fwhm * np.sqrt(4 * np.log(2) / np.pi) *
            np.exp(-4 * np.log(2) * ((x - cen) / fwhm)**2))


def lorentzian(x, area, cen, fwhm):
    return (2 * area / (np.pi * fwhm)) / (1 + 4 * ((x - cen) / fwhm)**2)


def pvoigt(x, area, cen, fwhm, frac):
    return frac * lorentzian(x, area, cen, fwhm) + (1 - frac) * gaussian(x, area, cen, fwhm)


# as a single function
def pvoigt_assym(x, area, cen, fwhm0, assym, frac):
    fwhm_x = (2 * fwhm0) / (1 + np.exp(assym * (x - cen)))
    return (frac * (2 * area / (np.pi * fwhm_x)) /
            (1 + 4 * ((x - cen) / fwhm_x)**2) +
            (1 - frac) * area / fwhm_x * np.sqrt(4 * np.log(2) / np.pi) *
            np.exp(-4 * np.log(2) * ((x - cen) / fwhm_x)**2))
