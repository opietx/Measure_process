import qcodes as qc
import numpy as np
from measure_process.fits.plots import plotter


from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import matplotlib.pyplot as plt
def Decay_formula(pars, x, data=None):
    amp = pars['amp'].value
    off = pars['off'].value
    T1 = pars['T1'].value
    
    model = amp*np.exp(-x/T1) + off

    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)
    
def fit_decay_T1(x,y, title='',plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=0.4, max=1.0, min=0.0)
    fit_params.add('off', value=0.7, max=1.0, min=0.0)
    fit_params.add('T1', value=x[int(len(x)/2)])
    
    out = minimize(Decay_formula, fit_params, args=(x,), kws={'data': y})
    
    
    if plot == True:
        plotter(x, y, out.params, title)
    print(f'T1 = {out.params["T1"].value} s')
    return out.params    
        
    
    
    

