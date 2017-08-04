import os
import re
import pickle
import pandas as pd


def elements():
    """
    Loads a DataFrame of all elements and isotopes.

    Returns
    -------
    pandas DataFrame with columns (element, atomic_number, isotope, atomic_weight, percent)
    """
    return pd.read_pickle(os.path.dirname(__file__) + '/periodic_table/elements.pkl')


def periodic_table():
    with open(os.path.dirname(__file__) + '/periodic_table/periodic_table.pkl', 'rb') as f:
        return pickle.load(f)


def calc_M(molecule):
    """
    Returns molecular mass of molecule.

    Where molecule is in standard chemical notation,
    e.g. 'CO2', 'HCO3' or B(OH)4

    """

    # load periodic table
    els = elements()
    els.set_index('element', inplace=True)

    # define regexs
    parens = re.compile('\(([A-z0-9]+)\)([0-9]+)?')
    stoich = re.compile('([A-Z][a-z]?)([0-9]+)?')

    ps = parens.findall(molecule)  # find subgroups in parentheses
    rem = parens.sub('', molecule)  # get remainder

    m = 0
    # deal with sub-groups
    if len(ps) > 0:
        for sub, ns in ps:
            ms = 0
            for e, n in stoich.findall(sub):
                me = (els.loc[e, 'atomic_weight'] *
                      els.loc[e, 'percent'] / 100).sum()
                if n == '':
                    n = 1
                else:
                    n = int(n)
                ms += me * n
            if ns == '':
                ns = 1
            else:
                ns = int(ns)
            m += ms * ns
    # deal with remainder
    for e, n in stoich.findall(rem):
        me = (els.loc[e, 'atomic_weight'] *
              els.loc[e, 'percent'] / 100).sum()
        if n == '':
            n = 1
        else:
            n = int(n)
        m += me * n
    return m


if __name__ == '__main__':
    print()
    print(calc_M('B(OH)3)'))
