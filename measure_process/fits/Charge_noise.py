
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


'''
Figure out units for reals
'''


from lmfit import Parameters, Minimizer, conf_interval, printfuncs


def one_over_f_noise(param, x, data=None):
    P_noise = param['P_noise']
    alpha = param['alpha']
    f0 = param['f0']
    B = param['B']
    
    # model = P_noise*x**(-alpha) + B/((x/f0)**2+1)
    model = P_noise*x**(-alpha)

    if data is None: 
        return model
    return np.nan_to_num((model - data)*1e6)
    


def fit_PSD(x,y,name='', plot=False):
    
    idx1 = ((x > 0.7) & (x <= 15)) #range
    
    #remove stable peaks
    idx2= ~((x>4.12) & (x<4.27)) 
    idx3= ~((x>5.51) & (x<5.8))
    idx4= ~((x>8.17) & (x<8.8))
    idx5= ~((x>11) & (x<11.46))
    idx6= ~((x>13.83) & (x<14))
    idx= np.where(idx1 & idx2 & idx3 & idx4 & idx5 & idx6) 
    
    
    x_fit, y_fit = x[idx], y[idx]
    # print(type(x_fit))
    
    
    
    fit_params = Parameters()
    fit_params.add('P_noise', value=1e-13)
    fit_params.add('alpha', value=1, min = 0.5, max=3)
    fit_params.add('B', value=1)
    fit_params.add('f0', value=1, min = 0.5, max=5, vary=False)
    
    
    mini = Minimizer(one_over_f_noise, fit_params, fcn_args=(x_fit, y_fit))
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
         fit = one_over_f_noise(out.params, x_fit)    
         
         plt.figure()
         plt.loglog(x_fit, y_fit)
         plt.plot(x_fit,fit,linewidth=3, alpha=0.8)
         
         plt.title(f'uuid={name} \n alpha = {out.params["alpha"].value} ')
         plt.ylabel('PSD')
         plt.xlabel('f [Hz')
         plt.show()
    # print(f'1HZ noise = {np.sqrt(10**out.params["P_noise"])}')
    print(f'alpha = {out.params["alpha"].value}')
    return out.params, errors
         

    
    
