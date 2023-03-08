import pandas as pd
import numpy as np
from PyAstronomy import pyasl

"""
load f1 table from NIST 
todo: formatting
"""

for i in range(1,118):
    AN = pyasl.AtomicNo()   
    atomic_symbol = AN.getElSymbol(atn=i)
    html_table = pd.read_html(f'https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z={i}&Formula=&gtype=0&range=U&lower=&upper=&density=&frames=no&htmltable=1')
    df = pd.DataFrame(np.array(html_table).reshape(-1,2), columns=['keV', 'f1'])
    df.to_csv(f'./f1/{i}_{atomic_symbol}_f1.csv')