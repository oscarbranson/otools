from .phreeq import input_str, run_phreeqc
import uncertainties as un
from scipy import stats

# Monte Carlo functions

class dummy_str:
    def __init__(self, string):
        self.string = string
    
    def rvs(self):
        return self.string

class dummy_numeric:
    def __init__(self, num):
        self.num = num
    
    def rvs(self):
        return self.num

def mc_input_dicts(input_dict, N, outputs=None):
    """
    Generates phreeqc input dicts for Monte-Carlo uncertainties
    
    Parameters
    ----------
    input_dict : dict 
        Containing `{entry: value}` pairs. The `value` may be either:
            - string, float or int : the value will be used for all iterations.
            - tuple : the first value will be taken as the mean, the second as 
              the standard deviation
            - uncertainties.core.Variable object : the nominal_value and std_dev
              will be used.
            - A scipy `rv_frozen` distribution object
            - a custom object that contains a `.rvs()` method to generate a draw
              for each iteration.
        In practice, all of the former are converted into objects with a 'rvs()'
        method before performing the iteration.
    N : int
        The number of monte-carlo iterations to generate.
    outputs : str or list
        a complete phreeqc output string, or a list of lines of an output string.

    Returns
    -------
    generator : where each iteration yields a new random dict drawn from the inputs.
    """
    
    # construct input dict of objects with .rvs() methods
    dists = {}
    for k, v in input_dict.items():
        if isinstance(v, str):
            dists[k] = dummy_str(v)
        elif isinstance(v, (float, int)):
            dists[k] = dummy_numeric(v)
        elif isinstance(v, un.core.Variable):
            dists[k] = stats.norm(v.nominal_value, v.std_dev)
        elif isinstance(v, tuple):
            dists[k] = stats.norm(v[0], v[1])
        elif hasattr(v, 'rvs'):
            dists[k] = v
        else:
            raise ValueError('Entry for {k} is invalid. See function doc for valid entry types.')
    
    for i in range(N):
        yield {k: v.rvs() for k, v in dists.items()}


def mc_input_str(input_dict, N, outputs=None):
    """
    Generates phreeqc input string for Monte-Carlo uncertainties
    
    Parameters
    ----------
    input_dict : dict 
        Containing `{entry: value}` pairs. The `value` may be either:
            - string, float or int : the value will be used for all iterations.
            - tuple : the first value will be taken as the mean, the second as 
              the standard deviation
            - uncertainties.core.Variable object : the nominal_value and std_dev
              will be used.
            - A scipy `rv_frozen` distribution object
            - a custom object that contains a `.rvs()` method to generate a draw
              for each iteration.
        In practice, all of the former are converted into objects with a 'rvs()'
        method before performing the iteration.
    N : int
        The number of monte-carlo iterations to generate.
    outputs : str or list
        a complete phreeqc output string, or a list of lines of an output string.

    Returns
    -------
    generator : where each iteration yields a new random dict drawn from the inputs.
    """

    return input_str(mc_input_dicts(input_dict=input_dict, N=N, outputs=outputs))