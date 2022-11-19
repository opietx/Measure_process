import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid


ds_noMW = load_by_uuid(1668611081604974076)
ds_MW = load_by_uuid(1668611106055974076)


x = ds_noMW['read12'].x()[10:-10]

noMW = np.nan_to_num(ds_noMW['read12'].y())[10:-10]
MW = np.nan_to_num(ds_MW['read12'].y())[10:-10]

# plt.figure()
# plt.plot(x,noMW, label='no_MW')
# plt.plot(x,MW, label='MW')
plt.plot(x,np.abs(MW-noMW), label='Vis_1')
plt.plot(x[np.argmax(np.abs(MW-noMW))],np.abs(MW-noMW).max(), 'o')
print(f'VIS={np.abs(MW-noMW).max()*100}%, vP1= {x[np.argmax(np.abs(MW-noMW))]} mV')
plt.legend()
plt.show()



def fit_fast_PSB(ds_noMW, ds_MW, measurement='read12', plot=True, fig_num=999):
    from scipy.signal import savgol_filter
    
    ds_noMW = load_by_uuid(ds_noMW)
    ds_MW = load_by_uuid(ds_MW)


    x = ds_noMW[measurement].x()[10:-10]
    noMW = np.nan_to_num(ds_noMW[measurement].y())[10:-10]
    MW = np.nan_to_num(ds_MW[measurement].y())[10:-10]
    noMW = savgol_filter(noMW, 9, 3) # window size 51, polynomial order 3
    MW = savgol_filter(MW, 9, 3) # window size 51, polynomial order 3
    
    if plot:
        plt.figure(fig_num)
        # plt.plot(x,noMW, label='no_MW')
        # plt.plot(x,MW, label='MW')
        plt.plot(x,np.abs(MW-noMW), label='Vis')
        plt.plot(x[np.argmax(np.abs(MW-noMW))],np.abs(MW-noMW).max(), 'o')
        print(f'VIS={np.abs(MW-noMW).max()*100}%, vP1= {x[np.argmax(np.abs(MW-noMW))]} mV')
    return x[np.argmax(np.abs(MW-noMW))]