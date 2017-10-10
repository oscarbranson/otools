"""
DePaolo's (2011) Surface Kinetic Model
"""
import numpy as np

def calc_Rb_m2(Rp, Rb0=6e-7):
    """
    Function for calculating Rb as a function of Rp.

    DePaolo's (2011) 'Model 2', described in Section 4.1.
    Varies Rb as a function of Rp**0.5 below a critical value (Rb0)

    Parameters
    ----------
    Rp : array-like
        Precipitation rate.
    Rb0 : float
        When Rp > Rb0, Rb = Rb0.
        When Rp < Rb0, Rb = Rb0 - fn(Rp**0.5)

    Returns
    -------
    Rb : array-like
    """
    Rb = np.full(Rp.shape, Rb0)
    ind = Rp < Rb0
    var = Rp[ind]**0.5
    var /= var.max()
    var *= Rb0
    Rb[ind] = var
    return Rb


def SKM(Rp, Kf, Keq, Rb=6e-7, mode=1):
    """
    DePaolo's (2011) Surface Kinetic Model
    
    Parameters
    ----------
    Kf : array-like
        The kinetic ('forward') fractionation factor.
        Kp will asymptote towards this value at high
        Rp.
    Keq : array-like
        The equilibrium fractionation factor. Kp will
        asymptote towards this values at low Rp.
    Rp : array-like
        Precipitation rate.
    Rb : array-like
        Ion detachment rate.
    mode : int
        Which mode to run the model in:
        1 = constant Rb
        2 = variable Rb, modified as a function of
            Rp**0.5 below Rb. Described in section
            4.1 of DePaolo (2011).

    Returns
    -------
    Kp : array-like
        Partitioning / fractionation of element in precipitated mineral.
    """
    if mode == 2:
        Rb = calc_Rb_m2(Rp, Rb)
    return Kf / (1 + Rb * (Kf / Keq - 1) / (Rp + Rb))


def yKp_SKM(xKp, xKf, xKeq, yKf, yKeq):
    """
    Calculate Kp of element y as a function of Kp of element x,
    given fractionation factors.

    Parameters
    ----------
    xKp : array-like
        Kp of the independent element.
    xKf, xKeq : array-like
        Kf and Keq of the independent element.
    yKf, yKeq : array-like
        Kf and Keq of the dependent element.

    Returns
    -------
    """
    return yKf / (1 + (xKf - xKp) * (yKf / yKeq - 1) / (xKp * (xKf / xKeq - 1)))
