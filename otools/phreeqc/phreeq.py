"""
Functions for calculating the speciation of C and B in precipitation solutions.
"""

import os
import pandas as pd
import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod

import pkg_resources as pkgrs
import os

def get_database_path(database_name='pitzer'):
    database_dir = os.path.join(pkgrs.resource_filename('otools.phreeqc', 'resources'), 'database')
    return os.path.join(database_dir, database_name.replace('.dat', '') + '.dat')

def make_solution(inputs, n=1):
    inp = [f"SOLUTION {int(n):d}"]
    for k, v in inputs.items():
        if isinstance(v, str):
            inp.append(f'    {k:20s}{v:s}')
        else:
            inp.append(f'    {k:20s}{float(v):.8e}')
    return '\n'.join(inp) + '\n'

def input_str(inputs, outputs=None):
    """
    Generate an input for calculating PHREEQC solutions.

    Parameters
    ----------
    inputs : dict, or list of dicts
        Where each key is a valid PHREEQC input key, and each 
        value is its value.
        
        If a dict of dicts, multiple solutions are specified for 
        calculation with the names of the input dicts.
    
    outputs : array-like
        A list of output lines.
    """
    # inputs
    solutions = []
    if isinstance(inputs, list):
        n = 0
        for v in inputs:
            solutions.append(make_solution(v, n))
            n += 1
    else:
        solutions.append(make_solution(inputs), 1)

    # outputs
    output = ['SELECTED_OUTPUT']
    if outputs is None:
        output.append("""    -pH
    -temperature
    -alkalinity
    -ionic_strength
    -totals Cl Na Mg K B Ca C S(6)
    -m OH- H+
    # molality outputs
    -m B(OH)4- B(OH)3 CaB(OH)4+ MgB(OH)4+ NaB(OH)4 B3O3(OH)4- B4O5(OH)4-2  # boron
    -m HCO3- CO3-2 CO2  # carbon
    -m SO4-2 HSO4-  # S
    # activities
    -a OH- H+
    -a B(OH)4- B(OH)3 CaB(OH)4+ MgB(OH)4+ NaB(OH)4 B3O3(OH)4- B4O5(OH)4-2  # boron
    -si Calcite Aragonite
    """)
    else:
        output += outputs
    
    return '\n'.join(solutions) + '\n' + '\n'.join(output) + '\nEND'

# def input_str(inputs, outputs):
#     """
#     Generate an input for calculating PHREEQC solutions.

#     Parameters
#     ----------
#     inputs : dict, or dict of dicts
#         Where each key is a valid PHREEQC input key, and each 
#         value is its value.
        
#         If a dict of dicts, multiple solutions are specified for 
#         calculation with the names of the input dicts.
    
#     outputs : dict
#         The output lines for phreeqc
#     """
#     pre = "SOLUTION 1\n"
    
#     inp = []
#     for k, v in inputs.items():
#         if isinstance(v, str):
#             inp.append('    {:20s}{:s}'.format(k, v))
#         else:
#             inp.append('    {:20s}{:.8e}'.format(k, v))
        
#     post = """SELECTED_OUTPUT
#     -pH
#     -temperature
#     -alkalinity
#     -ionic_strength
#     -totals Cl Na Mg K B Ca C S(6)
#     -m OH- H+
#     # minteq4 outputs
#     -m H3BO3 H2BO3- NaH2BO3 CaH2BO3+ H5(BO3)2- H8(BO3)3-  # boron
#     -m HCO3- NaHCO3 NaCO3- CO3-2 H2CO3 CaCO3 CaHCO3+  # carbon
#     # pitzer outputs
#     -m B(OH)4- B(OH)3 CaB(OH)4+ B3O3(OH)4- B4O5(OH)4-2  # boron
#     -m HCO3- CO3-2 CO2  # carbon
#     -si Calcite Aragonite
# END
#     """
    
#     return pre + '\n'.join(inp) + '\n' + post

def run_phreeqc(input_string, database=None, phreeq_path='/usr/local/lib/libiphreeqc.so', output_file=False):
    """
    Run input string in phreeqc with specified database.

    Parameters
    ----------
    input_string : str
        Valid phreeqc input string with SELECTED_OUTPUT.
    database : str
        Name of an included database to use (e.g. 'pitzer'), or 
        a complete path to a different phreeqc database (e.g. './path/to/pitzer.dat')
    phreeq_path : str
        Path to iphreeqc shared library. Defaults to '/usr/local/lib/libiphreeqc.so',
        which should work for standard installs on Linux machines

    Returns
    -------
    pandas.Series of calculated species
    """
    if database is None:
        print('No database specified  :  using pitzer')
        database = get_database_path()
    elif not os.path.exists(database):
        database = os.path.join(get_database_path(database))

    if not os.path.exists(database):
        raise ValueError(f"Can't phreeqc database: {database}\n   Please check that it exists.")

    phreeqc = phreeqc_mod.IPhreeqc(phreeq_path)
    phreeqc.load_database(database)
    if output_file:
        phreeqc.set_output_file_on()
    phreeqc.run_string(input_string)
    out = phreeqc.get_selected_output_array()
    return pd.DataFrame(out[1:], columns=out[0])

# def input_str_CB(temp=25, pH=8.1, Na=0, Cl=0, K=0, B=0, Ca=0, DIC=0, Mg=0, SO4=0, units='mol/L'):
#     """
#     Generate phreeqc input string for calculating C and B chemistry of solution.
#     """
#     template = """SOLUTION 1
#         temp      {temp:.3f}
#         pH        {pH:.3f}
#         pe        4
#         redox     pe
#         units     {units:s}
#         density   1
#         Cl        {Cl:.6f}
#         Na        {Na:.6f}
#         Mg        {Mg:.6f}
#         B         {B:.6f}
#         Ca        {Ca:.6f}
#         C         {DIC:.6f}
#         K         {K:.6f}
#         S(6)      {SO4:.6f}
#         -water    1 # kg

#     SELECTED_OUTPUT
#         -pH
#         -temperature
#         -alkalinity
#         -ionic_strength
#         -totals Cl Na Mg K B Ca C S(6)
#         -m OH- H+
#         # minteq4 outputs
#         -m H3BO3 H2BO3- NaH2BO3 CaH2BO3+ H5(BO3)2- H8(BO3)3-  # boron
#         -m HCO3- NaHCO3 NaCO3- CO3-2 H2CO3 CaCO3 CaHCO3+  # carbon
#         # pitzer outputs
#         -m B(OH)4- B(OH)3 CaB(OH)4+ B3O3(OH)4- B4O5(OH)4-2  # boron
#         -m HCO3- CO3-2 CO2  # carbon
#         -si Calcite Aragonite
#     END
#     """

#     return template.format(temp=temp,
#                            pH=pH,
#                            Na=Na,
#                            Cl=Cl,
#                            B=B,
#                            Ca=Ca,
#                            DIC=DIC,
#                            Mg=Mg,
#                            SO4=SO4,
#                            units=units,
#                            K=K)

# # function to calculate solution C and B speciation
# def calc_cb(temp=25, pH=8.1, Na=0, Cl=0, K=0, B=0, Ca=0, DIC=0, Mg=0, SO4=0, database='pitzer', phreeqc_path='/usr/local/lib/libiphreeqc.so', summ=True):
#     """
#     Calculate carbon and boron chemistry of solution.

#     Parameters
#     ----------
#     temp, pH, Na, Cl, B, Ca, DIC, Mg : float
#         Solution characteristics. All concentrations in mol/L.
#     database : str
#         Name or path of phreeqc database to use.
#     summ : boolean
#         If true, returns summary data, with column names tweaked to accommodate
#         different B species returned by pitzer / minteq.v4 databases.

#     Returns
#     -------
#     pandas.Series of results.
#     """
#     # create input string
#     inp = input_str(temp, pH, Na, Cl, K, B, Ca, DIC, Mg, SO4)
#     # run phreeqc
#     dat = run_phreeqc(inp, database, phreeqc_path)

#     if summ:
#         # return C and B chemistry (simple ions only)
#         out = pd.Series(index=['C', 'CO2', 'HCO3', 'CO3', 'B', 'BOH3', 'BOH4'])
#         # carbon species
#         out['C'] = dat['C(mol/kgw)']
#         if dat['m_CO2(mol/kgw)'] != 0:
#             out['CO2'] = dat['m_CO2(mol/kgw)']
#         else:
#             out['CO2'] = dat['m_H2CO3(mol/kgw)']

#         out['HCO3'] = dat['m_HCO3-(mol/kgw)']
#         out['CO3'] = dat['m_CO3-2(mol/kgw)']
#         # boron species
#         out['B'] = dat['B(mol/kgw)']
#         if dat['m_B(OH)3(mol/kgw)'] != 0:
#             out['BOH3'] = dat['m_B(OH)3(mol/kgw)']
#         else:
#             out['BOH3'] = dat['m_H3BO3(mol/kgw)']
#         if dat['m_B(OH)4-(mol/kgw)'] != 0:
#             out['BOH4'] = dat['m_B(OH)4-(mol/kgw)'] + dat['m_CaB(OH)4+(mol/kgw)']
#             out['BOH4_free'] = dat['m_B(OH)4-(mol/kgw)']
#         else:
#             out['BOH4'] = dat['m_H2BO3-(mol/kgw)'] + dat['m_NaH2BO3(mol/kgw)']  + dat['m_CaH2BO3+(mol/kgw)']
#             out['BOH4_free'] = dat['m_H2BO3-(mol/kgw)']

#         # auxillary data
#         out['pH'] = dat['pH']
#         out['temp'] = dat['temp(C)']
#         out['alk'] = dat['Alk(eq/kgw)']
#         out['SIc'] = dat['si_Calcite']
#         out['SIa'] = dat['si_Aragonite']
#         out['ion_str'] = dat['mu']

#         out['Ca'] = dat['Ca(mol/kgw)']
#         out['Na'] = dat['Na(mol/kgw)']
#         out['Cl'] = dat['Cl(mol/kgw)']
#         out['Mg'] = dat['Mg(mol/kgw)']
#         out['K'] = dat['K(mol/kgw)']
#         out['SO4'] = dat['S(6)(mol/kgw)']

#         return out
#     else:
#         return dat


# def calc_cb_rows(df, dbase='pitzer'):
#     """
#     Calculate solution conditions for each row of solution data.

#     Parameters
#     ----------
#     df : pandas.DataFrame
#         Each row must contain ['Temp (°C)', 'pH (NBS)', '[Na M)',
#         '[Cl M)', '[Ca M)', '[B M)', '[DIC M)', '[Mg M)]
    
#     Returns
#     -------
#     pandas.Dataframe with same index as input, with calculated C and B chemistry.
#     """
#     # determine column names
#     r = df.iloc[0,:]
#     cols = calc_cb(temp=r['Temp (°C)'],
#                    pH=r['pH (NBS)'],
#                    Na=r['[Na] (M)'], 
#                    Cl=r['[Cl] (M)'], 
#                    Ca=r['[Ca] (M)'], 
#                    B=r['[B] (M)'], 
#                    DIC=r['[DIC] (M)'], 
#                    Mg=r['[Mg] (M)'], dbase=dbase).index

#     # create empty output dataframe
#     out = pd.DataFrame(index=df.index, columns=cols)
#     out.sort_index(1, inplace=True)
#     out.sort_index(0, inplace=True)
    
#     for i, r in df.iterrows():
#         out.loc[i, :] = calc_cb(temp=r['Temp (°C)'],
#                                 pH=r['pH (NBS)'],
#                                 Na=r['[Na] (M)'], 
#                                 Cl=r['[Cl] (M)'], 
#                                 Ca=r['[Ca] (M)'], 
#                                 B=r['[B] (M)'], 
#                                 DIC=r['[DIC] (M)'], 
#                                 Mg=r['[Mg] (M)'], dbase=dbase)
#     return out
