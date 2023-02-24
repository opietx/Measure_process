import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid
from scipy import integrate
import scipy

from numba import jit


#%% 
plot_verbose = True
uuid = 1659132261101974087      
ds = load_by_uuid(uuid)

raw_init = (ds.m1_1()).ravel()
raw_read = (ds.m1_2()).ravel()
raw_total = np.concatenate((raw_init,raw_read))
# raw_total = raw_read


histogram  = np.histogram(raw_total, bins = int(np.sqrt(len(raw_total))/10), density=True)
hist=histogram[0]
bins=histogram[1]
bins = (bins[:-1] + bins[1:]) / 2

X,Y = bins, hist
plt.figure()
plt.hist(raw_total, bins = int(np.sqrt(len(raw_total))/10), alpha=0.4, label='data', density=True)
plt.show()
#%% initial fit 
from scipy.optimize import curve_fit

def gaussian(x,mu,sig):
    return (1/(np.sqrt(2*np.pi)*sig))*np.exp(-(x-mu)**2 /(2*sig**2))

def double_gaussian(x, prob_S, mu_1, sig_1, mu_2,sig_2):
    return prob_S*gaussian(x, mu_1, sig_1) + (1-prob_S)*gaussian(x, mu_2,sig_2)

# for SD2 triplets are on the right, singlets on the left (T=-80, S=-140)
X_center = abs((X.max()-X.min())/2) + X.min()
initial_params = [0.4,X_center-10, 8, X_center+10,9] # this is meant to be [mu_S,sig_S,mu_T,sig_T]
bound = ([0.1,X.min(),2, X.min(), 2],[0.8,X.max() , 10, X.max(), 10])
raw_fit, pcov = curve_fit(double_gaussian, X, Y, p0=initial_params, bounds = bound)

if plot_verbose:
    plt.plot(X,raw_fit[0]*gaussian(X,*raw_fit[1:3]), label = 'raw_fit S')
    plt.plot(X,(1-raw_fit[0])*gaussian(X,*raw_fit[3:]), label = 'raw_fit T')
plt.legend()

#%% fit with decays
@jit
def prob_decay_r(n,n_T):
    return np.exp(-n/n_T)
@jit
def integrand_T(s,x,mu_S,sig_S,mu_T,sig_T, n_T, n_r):
    return gaussian(x,mu_T,sig_T)*(prob_decay_r(n_r,n_T)/n_r)+ ((prob_decay_r(s,n_T)/n_T)*
            gaussian(x, (((s/n_r)*mu_T) +(((n_r-s)/n_r)*mu_S)), np.sqrt((((s/n_r)*sig_T)**2 +(((n_r-s)/n_r)*sig_S)**2))))


@jit
def integral(x,mu_S,sig_S,mu_T,sig_T, n_T,n_r):
    a, b= 0, n_r
    return scipy.integrate.quad(integrand_T,a,b, args=(x,mu_S,sig_S,mu_T,sig_T, n_T, n_r))[0]

@jit
def C_T(x,mu_S,sig_S,mu_T,sig_T, n_T, n_r):
    cc =np.empty((len(x)))
    for ii, x_s in  enumerate(x):
        rr = integral(x_s,mu_S,sig_S,mu_T,sig_T, n_T,n_r)
        cc[ii] = rr
    return cc


    


#%% fit histogram
T1=80e-6
sampling_rate = 1e9
t=40e-6
n_T = T1/(1/sampling_rate)
n_samples = round(sampling_rate*t)
n_r = n_samples

def model_final_integral(x,prob_S,mu_S,sig_S,mu_T,sig_T, T1):
    n_T = T1/(1/sampling_rate)
    return prob_S*gaussian(x,mu_S,sig_S) + (1-prob_S)*C_T(x,mu_S,sig_S,mu_T,sig_T, n_T, n_r)


init_vals = np.concatenate((raw_fit,[T1]))
bound = ([0.2,X.min(),3, X.min(), 3,5e-6],[0.8,X.max() , 10, X.max(), 10,1e-3])
popt, pcov = scipy.optimize.curve_fit(model_final_integral, X, Y, p0=init_vals, bounds = bound)
plt.plot(X,model_final_integral(X,*popt), label = 'Final Fit', alpha=0.7, linewidth=4)
plt.legend()
plt.show()
perr = np.sqrt(np.diag(pcov))
if plot_verbose:
    n_T = popt[5]/(1/sampling_rate)
    fit_triplet_decay = (1-popt[0])*C_T(X,*popt[1:5],n_T, n_samples)
    fit_singlet  = popt[0]*gaussian(X,*popt[1:3])
    plt.plot(X,fit_singlet, label = 'fit_S')
    plt.plot(X,fit_triplet_decay, label = 'fit_T')
    plt.legend()

    for ii in range(len(popt)):
        print(f'{round(popt[ii],4)} +- {round(perr[ii],4)}')
print(f'T1 time = {popt[-1]*1e6}us')

#%%

#%% Fidelities  functions
@jit
def F_T(x,mu_S,sig_S,mu_T,sig_T, n_T,n_r):
    cc =np.empty((len(x)))
    for ii,x_s in enumerate(x):
        rr = integral2(x_s,mu_S,sig_S,mu_T,sig_T, n_T,n_r)
        cc[ii] = rr
    return cc
@jit
def integral2(x,mu_S,sig_S,mu_T,sig_T, n_T,n_r):
    a1, b1= 0, n_r
    a2, b2= 0, x
    return scipy.integrate.dblquad(integrand_T,a2,b2,a1,b1, args=(mu_S,sig_S,mu_T,sig_T, n_T, n_r))[0]
@jit
def F_gaussian(x,mu,sig):
    return 0.5*(1+scipy.special.erf((x-mu)/(sig*np.sqrt(2))))


#%%
n_T = popt[5]/(1/sampling_rate)
fid_T = F_T(X,*popt[1:5],n_T,n_samples)+1
fid_S = F_gaussian(X, *popt[1:3])

V_E = abs(fid_S-fid_T)
v_optimized = X[np.argmax(V_E)]

plt.figure()
plt.plot(X,fid_T, '-',label=r'$F_{|T\rangle}^{decay}$')
plt.plot(X,fid_S, '-',label=r'$F_{|S\rangle}$')
plt.plot(X,V_E, '-',alpha=0.5,linewidth=5,label=r'$V^{decay}$')
plt.vlines(v_optimized,0,max(V_E),colors=['green'], linestyles='dashed')
plt.hlines(max(V_E),v_optimized,X[0],colors=['green'], linestyles='dashed')
plt.xlabel('Current [A]')
plt.ylabel('Fidelities/visibility')

tit = r'$V_{max}^{decay}$'
tit =tit+ '={}%\n'.format(round(np.max(V_E)*100,2))
plt.title(tit)
plt.legend()
plt.show()



