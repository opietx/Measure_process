import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from measure_process.Init_file import data_dump, query_database
from core_tools.data.ds.data_set import load_by_uuid

from measure_process.fits.T2_star import  Ramsey_formula, fit_Ramsey
from measure_process.fits.Readout_point import fit_fast_PSB
from measure_process.fits.Peaks import  fit_Gaussian

def find_closest_idx(arr, val):
       idx = np.abs(arr - val).argmin()
       return idx
   
def process_data(ds_avg, mask):
    T2 = []
    omega = []
    flag=False
    for ii,data in tqdm(enumerate(ds_avg), total=(ds_avg.shape[0]-len(mask)), desc='Fitting'):
        if ii in mask:
            pass
        else:
            try:
                params, errors = fit_Ramsey(1e-9*data.x()[:], data.y()[:], uuid, plot=plot)
                T2.append(params['T2'].value)
                params, errors = fit_Ramsey(1e-9*data.x()[0:70], data.y()[0:70], uuid, plot=plot)
                omega.append(params['omega'].value)
            except:
                pass
    return np.array(T2), np.array(omega)


def nan_mask(ds_avg):   
    '''
    selects traces in which all points are nan, aka have not been measured

    Parameters
    ----------
    ds_avg : TYPE
        DESCRIPTION.

    Returns
    -------
    mask : TYPE
        indices where we have nan everywhere.

    '''
    mask = []
    for ii,data in enumerate(ds_avg):
        if np.isnan(data.y()).all():
            mask.append(ii)
    return mask

def vals_mask(ds_avg,thr):    
    '''
    Values that i want to throw away because the charge jump.
    If thr=0 then we keep all the data
    Parameters
    ----------
    ds_avg : qcodes_measurement
        
    thr : double
        threshold where to make the cut, we keep data below it

    Returns
    -------
    mask : array
        contains the indices of data to throw away

    '''
    mask = []
    for ii,data in enumerate(ds_avg):
        if ds_avg.x()[ii]>thr:
            mask.append(ii)
    if thr == 0:
        mask=[]
    return mask

#%% obtain the readout traces
datasets = query_database('2022-11-16 10:50:30','2022-11-16 16:20:30', name = 'MW' )
for ii in tqdm(range(int(len(datasets['uuid'])/2 -1 ))):
    ii *= 2
    noMW = datasets['uuid'][ii]
    MW = datasets['uuid'][ii+1]
    ds_noMW = load_by_uuid(noMW)
    ds_MW = load_by_uuid(MW)
    fit_fast_PSB(noMW,MW)
    


#%% 1 dataset averaged all
qubit = 1
label='vP1 = 21.3mV'
plot = False 
uuids = [1668609555987974076]

#%%
T2s, Omegas = np.array([]),np.array([])
for uuid in uuids:
    ds = load_by_uuid(uuid)
    ds_avg = ds[f'read{qubit}']
       
    mask1 = nan_mask(ds_avg)
    mask2 = vals_mask(ds_avg, 0)#introduce X axis value for thresholding when the charge jump is. keep data below it
    mask = np.unique(np.concatenate((mask1,mask2),0))

    T2, Omega = process_data(ds_avg, mask)
    T2s = np.append(T2s,T2)
    Omegas = np.append(Omegas,Omega)

#%% histogram T2*
name = f'{label}, N={len(T2s)}'
counts, bins = np.histogram(T2s, density=True, range=(2e-6,3.5e-6), bins= int(3*np.sqrt(len(T2s))))
params, errors = fit_Gaussian(bins,counts, title=name, plot=True)
#%% 
param = Omegas


dt = (ds_avg.x()[1]-ds_avg.x()[0])
time = np.arange(0,dt*len(param), dt)
freq = np.fft.rfftfreq(len(time), d=dt)
fourier = np.abs(np.fft.rfft(param))

fig, ax = plt.subplots(2,num=99)
ax[0].set_ylabel('Qubit freq detuning [Hz]')
ax[0].set_xlabel('time [s]')
ax[0].plot(time, param)

ax[1].loglog(freq,fourier,label=label)
ax[1].set_xlabel('freq [Hz]')
ax[1].set_ylabel('PSD [Hz**2/Hz]')
handles, labels = ax[1].get_legend_handles_labels()
fig.legend(handles, labels)
fig.show()

#%% plot Omega with datasets and charge jumps 
plt.figure(900)
x_last = 0
qubit=1
plot = False 
datasets = query_database('2022-11-19 8:53:30','2022-11-19 9:08:30', name = 'Ramsey_time_evo' )
thresholds = [160,0,400,80,0,500]#this are in minutes 
color = 'blue'
for ii,uuid in enumerate(datasets['uuid']):
    ds = load_by_uuid(int(uuid))
    ds_avg = ds[f'read{qubit}']
    
    thr = thresholds[ii]
    thr_idx = find_closest_idx(ds_avg.x(),thr)
    print(thr_idx)
    
    
       
    mask1 = nan_mask(ds_avg)
    mask2 = vals_mask(ds_avg, 0)#introduce X axis value for thresholding when the charge jump is. keep data below it
    mask = np.unique(np.concatenate((mask1,mask2),0))

    T2, Omega = process_data(ds_avg, mask)
    # print(Omega)
    dt = (ds_avg.x()[1]-ds_avg.x()[0])
    time = x_last+np.arange(0,dt*len(Omega), dt)
    if (thr == 0) or (thr>dt*len(Omega)):
        plt.plot(time,Omega, color = color)
    else:
        plt.plot(time[:thr_idx],Omega[:thr_idx], color = color)
        if color=='blue':
            color='green'
        else:
            color='blue'
        plt.plot(time[thr_idx+1:],Omega[thr_idx+1:], '--' ,color=color)
        plt.axvline(x_last+thr, linewidth = 1, alpha = 0.3, color='red')
        
    x_last = time[-1]
    plt.axvline(x_last, linewidth = 1, alpha = 0.3, color='gray')
    
# plt.legend()
plt.show()
    
    
