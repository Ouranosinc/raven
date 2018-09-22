from raven import models
import pandas as pd
import numpy as np



def test_gr4j():
    time = pd.date_range('2000-01-01', freq='D', periods=365*3)
    n = len(time)
    temp = pd.Series(-20*np.sin(np.arange(n)/365.25 * 2 * np.pi)+5, index=time)
    precip = pd.Series(((np.random(n)-.5)*10).clip(0), index=time)
    evap = np.clip(temp/200 + np.random(n)*10, 0)
    df = pd.DataFrame({'Temp':temp, 'Evap':evap, 'Prec':precip})

    return df
