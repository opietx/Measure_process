import qcodes as qc
import numpy as np
import matplotlib.pyplot as plt


from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import matplotlib.pyplot as plt
def Rabi_formula(pars, x, data=None):
    amp = pars['amp'].value
    off = pars['off'].value
    w1 = pars['freq'].value
    Dw = pars['detuning'].value
    
    
    model = amp*w1**2/(w1**2+Dw**2)*np.sin(x/2 * np.sqrt(w1**2+Dw**2))**2 +off

    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)
    
def fit_Rabi(x,y, title='',plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=1, max=1, min=0.999)
    fit_params.add('off', value=y.min(), max=1, min=0)
    fit_params.add('freq', value=15e6)
    fit_params.add('detuning', value=-6e6)
    
    out = minimize(Rabi_formula, fit_params, args=(x,), kws={'data': y})
    
    
    if plot == True:
        fit = Rabi_formula(out.params, x)    
        plt.figure()
        plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        
        plt.ylabel('fraction spin-up')
        plt.xlabel('Rabi time [s]')
        plt.show()
    print(f'detuning = {out.params["detuning"].value} ')
    return out.params    
        