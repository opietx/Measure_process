import qcodes as qc
import numpy as np
import matplotlib.pyplot as plt


from lmfit import Parameters, Minimizer, conf_interval, printfuncs

def Decay_formula(pars, x, data=None):
    A_uu = pars['A_uu'].value
    A_ud = pars['A_ud'].value
    A_du = pars['A_du'].value
    off = pars['off'].value
    
    T1_q1 = pars['T1_q1'].value
    T1_q2 = pars['T1_q2'].value
    
    model = (A_uu+A_ud)*np.exp(-x/T1_q1) + (A_uu+A_du)*np.exp(-x/T1_q2) - 2*A_uu*np.exp(-x/T1_q2)*np.exp(-x/T1_q1) + off

    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)
    
def fit_mixed_decay_T1(x,y, title='',plot=False):
    fit_params = Parameters()
    fit_params.add('A_uu', value=0.1, max=1.0, min=0)
    fit_params.add('A_ud', value=0.4, max=1.0, min=0)
    fit_params.add('A_du', value=0.4, max=1.0, min=0)

    fit_params.add('off', 0.1, max=1.0, min=-1.0)
    fit_params.add('T1_q1', value=x[int(len(x)/2)], min=0)
    fit_params.add('T1_q2', value=x[int(len(x)/2)], min=0)
    print(fit_params)
    mini = Minimizer(Decay_formula, fit_params, fcn_args=(x, y))
    out = mini.minimize()
    errors = {}
    # try:
    #     ci =  conf_interval(mini, out, sigmas=[2], trace=False)
    #     for item in out.params.keys():
    #         errors[item] = [out.params[item].value-ci[item][0][1],ci[item][2][1]-out.params[item].value]    
    # except:
    #     for item in out.params.keys():
    #         errors[item] = [0,0]    
    # printfuncs.report_ci(ci)
    
    
    
    if plot == True:
        fit = Decay_formula(out.params, x)    
        # plt.figure()
        plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        plt.title(f'uuid={title} \n T1_q1 = {out.params["T1_q1"].value} s \n T1_q2 = {out.params["T1_q2"].value} s')
        plt.ylabel('fraction spin-up')
        plt.xlabel('waiting time [s]')
        plt.show()
    # print(f'T1 = {out.params["T1"].value} s')
    return out.params, errors
        
    
    
#%%
qubit = 1
plot = True 
uuid = 1652890383297974076   
ds = load_by_uuid(uuid)
# data_plotter(ds)
ds_avg = ds[f'read{qubit}'].average('x')
# ds_avg_y = ds_avg.y()-ds_avg.y().min()
# ds_avg_y /= ds_avg_y.max()
#it is 1-y for T1 and only y for T1_PSB
params, errors = fit_mixed_decay_T1(1e-9*ds_avg.x(), 1-ds_avg.y(), uuid, plot=plot)
