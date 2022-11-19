import numpy as np
import matplotlib.pyplot as plt
from numba import jit


from lmfit import Parameters, Minimizer, conf_interval, printfuncs
def Ramsey_formula(pars, x, data=None):
    amp = pars['amp'].value
    off = pars['off'].value
    omega = pars['omega'].value
    phi = pars['phi'].value
    T2 = pars['T2'].value
    alpha = pars['alpha'].value

    model = amp*np.cos(2*np.pi*omega*x + phi)*np.exp(-(x/T2)**alpha) + off
    
    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)

# @jit # Set "nopython" mode for best performance, equivalent to @njit
def fit_Ramsey(x,y, title='', plot=False,freq_init=2e6):
    fit_params = Parameters()
    fit_params.add('amp', value=0.3, max=1.0, min=-1.0)
    fit_params.add('off', value=0.5, max=1.0, min=-1.0)
    fit_params.add('omega', value=freq_init)
    fit_params.add('phi', value=0)
    fit_params.add('T2', value=5e-6, min = 0)
    fit_params.add('alpha', value=2, min = 1, max = 3, vary=False)
    
    mini = Minimizer(Ramsey_formula, fit_params, fcn_args=(x, y))
    out = mini.minimize()
    errors = {}
    try:
        ci =  conf_interval(mini, out, sigmas=[2], trace=False)
        for item in out.params.keys():
            errors[item] = [out.params[item].value-ci[item][0][1],ci[item][2][1]-out.params[item].value]    
    except:
        for item in out.params.keys():
            errors[item] = [0,0]    
    # printfuncs.report_ci(ci)
    
    
    
    if plot == True:
        fit = Ramsey_formula(out.params, x)    
        
        # plt.figure()
        plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        # fit = Ramsey_formula(fit_params, x)    
        # plt.plot(x,fit,'-')
        plt.title(f'uuid={title} \n T2 = {out.params["T2"].value*1e6} us')
        plt.ylabel('fraction spin-up')
        plt.xlabel('waiting time [s]')
        plt.show()
    # print(f'T2 = {out.params["T2"].value*1e6} us')
    return out.params, errors
        

def fit_CPMG(x,y, title='', plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=0.4, max=1.0, min=-1.0)
    fit_params.add('off', value=0.4, max=1.0, min=-1.0)
    fit_params.add('omega', value=2e6)
    fit_params.add('phi', value=0)
    fit_params.add('T2', value=40e-6, min=1e-7) 
    
    mini = Minimizer(Ramsey_formula, fit_params, fcn_args=(x, y))
    out = mini.minimize()
    errors = {}
    try:
        ci =  conf_interval(mini, out, sigmas=[2], trace=False)
        for item in out.params.keys():
            errors[item] = [out.params[item].value-ci[item][0][1],ci[item][2][1]-out.params[item].value]    
    except:
        for item in out.params.keys():
            errors[item] = [0,0]    
    # printfuncs.report_ci(ci)
    
    
    
    if plot == True:
        # plt.figure()
        plt.plot(x,y,'o')
        
        fit = Ramsey_formula(out.params, x)    
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        # fit = Ramsey_formula(fit_params, x)    
        # plt.plot(x,fit,'-')
        plt.title(f'uuid={title} \n T2 = {out.params["T2"].value*1e6} us')
        plt.ylabel('fraction spin-up')
        plt.xlabel('waiting time [s]')
        plt.show()
    print(f'T2 = {out.params["T2"].value*1e6} us')
    return out.params, errors
        