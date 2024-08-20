import pandas as pd
import uncertainty_helpers as un
import uncertainties.unumpy as unp

def df_separate_uncertainties(df):
    """Create a copy of a dataframe where all columns containing uncertainties are split into two columns, one for the nominal value and one for the standard deviation.

    Parameters
    ----------
    df : pandas.DataFrame
        A pandas dataframe where some of the columns are uncertainties.unumpy.uarrays.
        
    Returns
    -------
    A copy of the dataframe with an additional level added to the column index identifying the 'mean' and 'std' values.
    """
    
    sdf = pd.DataFrame(index=df.index, columns=pd.MultiIndex(levels=[[],[],[]], codes=[[],[],[]], names=['Element', 'Unit', 'Type']))

    for c, d in df.items():
        if isinstance(d.iloc[0], un.core.AffineScalarFunc):
            sdf.loc[:, (*c, 'mean')] = unp.nominal_values(d)
            sdf.loc[:, (*c, 'std')] = unp.std_devs(d)
        else:
            sdf.loc[:, (*c, '')] = d
     
    return sdf
    