import importlib.resources as resources
import pandas as pd

package_path = resources.files('otools')

def load():
    f = package_path / 'seawater/seawater.csv'
    
    return pd.read_csv(f, comment='#').set_index('name')