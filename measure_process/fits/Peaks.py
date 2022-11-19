import numpy as np
import matplotlib.pyplot as plt


from lmfit import Parameters, Minimizer, conf_interval, printfuncs

def gaussian_formula(pars, x, data=None):
    sigma = pars['sigma'].value
    avg = pars['avg'].value
    amp = pars['amp'].value

    model = 1/np.sqrt(sigma**2 * 2*np.pi)*np.exp(-0.5*(x-avg)**2/sigma**2)*amp
    
    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)


def fit_Gaussian(x,y, title='', plot=False):
    fit_params = Parameters()
    fit_params.add('sigma', value=x[int(len(x)/2)]/5)
    fit_params.add('amp', value=1, min = 0, max = y.max())
    fit_params.add('avg', value=x[int(len(x)/2)], max = x[-1], min=x[0])
    x_centers = (x[:-1] + x[1:]) / 2
    
    mini = Minimizer(gaussian_formula, fit_params, fcn_args=(x_centers, y))
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
        fit = gaussian_formula(out.params, x_centers)    
        
        # plt.figure()
        plt.stairs(y, x, alpha = 0.5, label=f'{title}', fill=True)
        
        plt.plot(x_centers,fit,linewidth=5, alpha=0.5)
        # fit = gaussian_formula(fit_params, x_centers)    
        # plt.plot(x_centers,fit,'-')
        plt.ylabel('counts [#]')
        plt.xlabel('T2* [s]')
        plt.legend()
        plt.show()
        
    return out.params, errors