def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


#%% Plot SD and thr form the SD_calibration
ds = load_by_uuid(1655294163383974076)
plt.figure()
subset = ds.m1_1.average(2)
data = subset()
y = subset.y()

idx = np.argmax(data[0] - data[1])
SD1_winner = y[idx]
threshold = np.average(data[:,idx])

if True:
    plt.plot(y, data[0])
    plt.plot(y, data[1])
    plt.scatter(SD1_winner, threshold)
    
    plt.xlabel('SD volatage (mV)')
    plt.ylabel('Reflected signal (mV)')
    plt.show()

print(f'\n\n--SD1 calib complete--\n\n\tSetting SD1_P to : {SD1_winner},\n\tSetting threshold set to {threshold}')

#%%

ds = load_by_uuid(1655297149004974076)
data_plotter(ds)



data = ds['RAW read1 (SD1_IQ:1)'].z()
data = data[:].flatten()


plt.figure()
hist = plt.hist(data,bins = int(np.sqrt(len(data))/5), label='Readout signal')
plt.vlines(-127.5, 0, hist[0].max(), label='signal hreshold', alpha=0.5, linewidth=5)


plt.xlabel('Signal [mV]')
plt.ylabel('Counts')
plt.title('Readout time 10us')
plt.legend()
plt.show()

#%%

from scipy.signal import savgol_filter

ds = load_by_uuid(1655294877875974076)
# plt.figure()
data = ds['RAW init1 (SD1_IQ:0)'][550:650]
data = data.z()
data = data.flatten()



hist = plt.hist(data,bins = int(np.sqrt(len(data))/3))
# plt.vlines(-132, 0,np.max(hist[0]), alpha=0.5, linewidth=5)

plt.xlabel('Signal [mV]')
plt.ylabel('Counts')
plt.title('Readout time 70us')
# plt.legend()
plt.show()


#%%
from scipy.signal import savgol_filter

ds = load_by_uuid(1655294219136974076)
data_plotter(ds)
step = 60
maxs = []
plt.figure()
for ii in range(0,int(600/step)):
    data = ds['RAW init1 (SD1_IQ:0)'][int(ii*step):int((ii+1)*step)]
    data = data.z()
    data = data.flatten()
    
    
    hist = np.histogram(data,bins = int(np.sqrt(len(data))/3))
    yhat = hist[0]
    # yhat = savgol_filter(yhat, 21, 3)
    plt.plot(hist[1][:-1],yhat)
    
    x = hist[1][int(len(yhat)/2):]
    maxs.append(x[np.argmax(yhat[int(len(yhat)/2):])])


plt.xlabel('Signal [mV]')
plt.ylabel('Counts')
plt.title('Steps of 500ns')
# plt.legend()
plt.show()


temps = (np.linspace(int(step/2),int(600-step),int(600/step))*ds['RAW init1 (SD1_IQ:0)'].x()[1])
maxs = np.array(maxs)
plt.figure()
plt.plot(temps/1000,maxs,'o')
plt.xlabel('pulse_duration [us]')
plt.ylabel('max_second peak')
# plt.title('Steps of 500ns')



from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import matplotlib.pyplot as plt
def Decay_formula(pars, x, data=None):
    amp = pars['amp'].value
    off = pars['off'].value
    T1 = pars['T1'].value
    
    model = 1-amp*np.exp(-x/T1) - off

    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)
    
def fit_decay_T1(x,y, title='',plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=abs(y.max()-y.min()))
    fit_params.add('off', value=y.min())
    fit_params.add('T1', value=x[int(len(x)/2)])
    
    out = minimize(Decay_formula, fit_params, args=(x,), kws={'data': y})
    
    
    if plot == True:
        fit = Decay_formula(out.params, x)    
        # plt.figure()
        # plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        # plt.title(f'uuid={title} \n T1 = {out.params["T1"].value} s')
        # plt.ylabel('fraction spin-up')
        # plt.xlabel('waiting time [s]')
        # plt.show()
    print(f'T1 = {out.params["T1"].value} s')
    return out.params    

params = fit_decay_T1(temps/1000, maxs, 'test', plot=True)
