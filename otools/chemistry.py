"""
The periodic table, and all it's info! And functions for doing chemical things.
"""
import os
import re
import pickle
import pandas as pd


def elements(all_isotopes=True):
    """
    Loads a DataFrame of all elements and isotopes.

    Scraped from https://www.webelements.com/

    Returns
    -------
    pandas DataFrame with columns (element, atomic_number, isotope, atomic_weight, percent)
    """
    el = pd.read_pickle(os.path.dirname(__file__) + '/periodic_table/elements.pkl')
    if all_isotopes:
        return el
    else:
        def wmean(g):
            return (g.atomic_weight * g.percent).sum() / 100
        iel = el.groupby('element').apply(wmean)
        iel.name = 'atomic_weight'
        return iel


def periodic_table():
    """
    Loads dict containing all elements and associated metadata.

    Scraped from https://www.webelements.com/

    Returns
    -------
    dict
    """
    with open(os.path.dirname(__file__) + '/periodic_table/periodic_table.pkl', 'rb') as f:
        return pickle.load(f)

def decompose_molecule(molecule, n=1):
    """
    Returns the chemical constituents of the molecule, and their number.

    Parameters
    ----------
    molecule : str
        A molecule in standard chemical notation, 
        e.g. 'CO2', 'HCO3' or 'B(OH)4'.
    
    Returns
    -------
    All elements in molecule with their associated counts : dict
    """
    if isinstance(n, str):
        n = int(n)
    
    # define regexs
    parens = re.compile('\(([A-z0-9()]+)\)([0-9]+)?')
    stoich = re.compile('([A-Z][a-z]?)([0-9]+)?')

    ps = parens.findall(molecule)  # find subgroups in parentheses
    rem = parens.sub('', molecule)  # get remainder
    
    if len(ps) > 0:
        for s, ns in ps:
            comp = decompose_molecule(s, ns)
        for k, v in comp.items():
            comp[k] = v * n
    else:
        comp = {}
        
    for e, ns in stoich.findall(rem):
        if e not in comp:
            comp[e] = 0
        if ns == '':
            ns = 1 * n
        else:
            ns = int(ns) * n
        comp[e] += ns

    return comp

def calc_M(molecule):
    """
    Returns molecular weight of molecule.

    Parameters
    ----------
    molecule : str
        A molecule in standard chemical notation, 
        e.g. 'CO2', 'HCO3' or 'B(OH)4'.
    
    Returns
    -------
    Molecular weight of molecule : dict
    """
    # load periodic table
    els = elements(all_isotopes=False)

    comp = decompose_molecule(molecule)

    m = 0
    for k, v in comp.items():
        m += els[k] * v
    
    return m

def seawater(Sal=35., unit='mol/kg'):
    """
    Standard mean composition of seawater.

    From Dickson, Sabine and Christian (2007), Chapter 5, Table 3

    @book{dickson2007guide,
          title={Guide to best practices for ocean CO2 measurements.},
          author={Dickson, Andrew Gilmore and Sabine, Christopher L and Christian, James Robert},
          year={2007},
          publisher={North Pacific Marine Science Organization},
          howpublished="https://www.nodc.noaa.gov/ocads/oceans/Handbook_2007.html",
          ISBN="1-897176-07-4"}
        

    Parameters
    ----------
    Sal : float
        Salinity, default is 35
    unit : str
        Either 'mol/kg' or 'g/kg'.

    Returns
    -------
    Seawater composition in chosen units at specified salinity : dict
    """
    
    sw = {"Cl": 0.54586,
          "SO4": 0.02824,
          "Br": 0.00084,
          "F": 0.00007,
          "Na": 0.46906,
          "Mg": 0.05282,
          "Ca": 0.01028,
          "K": 0.01021,
          "Sr": 0.00009,
          "B": 0.00042}

    for s in sw.keys():
        sw[s] *= Sal / 35.
    
    if unit == 'g/kg':
        for k, v in sw.items():
            sw[k] = calc_M(k) * v
    
    return sw

if __name__ == '__main__':
    print()
    print(calc_M('B(OH)3)'))
