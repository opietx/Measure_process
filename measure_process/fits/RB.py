from lmfit import Parameters, Minimizer, conf_interval, printfuncs
import matplotlib.pyplot as plt
import numpy as np

def rb_formula(pars,x, data=None):
    A = pars['A']
    B = pars['B']
    F = pars['Fid']
    
    model = A*(2*F-1)**x + B
    if data is None:
        return model
    
    return np.nan_to_num(model - data)
    
def fit_RB(x,y, plot=False):
    fit_params = Parameters()
    fit_params.add('A', value=-0.3, max=0.5) # read out fidelity/2
    fit_params.add('B', value=0.4, min=0.1,max=0.9) # offset of the visibility
    fit_params.add('Fid', value=0.99)
    
    mini = Minimizer(rb_formula, fit_params, fcn_args=(x, y))
    out = mini.minimize()
    errors = {}
    try:
        ci =  conf_interval(mini, out, sigmas=[2], trace=False)
        for item in out.params.keys():
            errors[item] = [out.params[item].value-ci[item][0][1],ci[item][2][1]-out.params[item].value]    
    except:
        for item in out.params.keys():
            errors[item] = [0,0]    
    
    
    if plot == True:
        #guess
        # fit = rb_formula(fit_params, x)    
        # plt.plot(x,fit,'-')
        
        #fit
        fit = rb_formula(out.params, x)
        plt.semilogx(x,y,'o')
        plt.semilogx(x,fit)
        plt.show()
    # print(f'fidelity = {out.params["Fid"].value*100}%')
    print(f'single gate fidelity = {(1-(1-out.params["Fid"].value)/1.875)*100} +- {out.params["Fid"].stderr*100/1.875}%')
    
    return out.params, errors

