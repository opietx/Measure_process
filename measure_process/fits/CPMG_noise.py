import qcodes as qc
import numpy as np
import matplotlib.pyplot as plt


from lmfit import Parameters, Minimizer, conf_interval, printfuncs
def beta_formula(pars, x, data=None):
    amp = pars['amp'].value
    # off = pars['off'].value
    beta = pars['beta'].value

    model = amp*x**(beta/(1+beta)) 
    
    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)
    
def fit_Beta(x,y, title='', plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=1)
    # fit_params.add('off', value=0.3)
    fit_params.add('beta', value=1, max=3, min=-3)
    
    mini = Minimizer(beta_formula, fit_params, fcn_args=(x, y))
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
        fit = beta_formula(out.params, x)    
        
        # plt.figure()
        plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        # fit = Ramsey_formula(fit_params, x)    
        # plt.plot(x,fit,'-')
        plt.title(f'uuid={title} \n B = {out.params["beta"].value}')
        plt.ylabel('T2_CPMG [s]')
        plt.xlabel('N-pulses')
        plt.show()
    print(f'beta = {out.params["beta"].value}')
    return out.params, errors
        