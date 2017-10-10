"""
Convenience functions for handling uncertainties.
"""

import uncertainties as un
import uncertainties.unumpy as up

def err(x):
    return up.std_devs(x)


def nom(x):
    return up.nominal_values(x)