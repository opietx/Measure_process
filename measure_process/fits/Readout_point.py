import matplotlib.pyplot as plt
import numpy as np


def fit_fast_PSB(ds_noMW, ds_MW, measurement='read12', plot=True, fig_num=999):
    '''
    Plots the difference between MW_on and MW_off and finds the spot with highest
    visibility

    Parameters
    ----------
    ds_noMW : TYPE
        dataset, loaded
    ds_MW : TYPE
        dataset, loaded
    measurement : TYPE, optional
        measurement used to read from . The default is 'read12'.
    plot : TYPE, optional
        show plot. The default is True.
    fig_num : TYPE, optional
        so that eveyrthing appears on the same plot. The default is 999.

    Returns
    -------
    readout_point : double
        mV of the readout point with max visibility.
    vis : double
        visibility of that point

    '''
    from scipy.signal import savgol_filter
    

    x = ds_noMW[measurement].x()[10:-10]
    noMW = np.nan_to_num(ds_noMW[measurement].y())[10:-10]
    MW = np.nan_to_num(ds_MW[measurement].y())[10:-10]
    noMW = savgol_filter(noMW, 9, 3) # window size 51, polynomial order 3
    MW = savgol_filter(MW, 9, 3) # window size 51, polynomial order 3
    
    readout_point, vis = x[np.argmax(np.abs(MW-noMW))], np.abs(MW-noMW).max()
    if plot:
        plt.figure(fig_num)
        # plt.plot(x,noMW, label='no_MW')
        # plt.plot(x,MW, label='MW')
        plt.plot(x,np.abs(MW-noMW), alpha=0.5, linewidth = 3,label='Vis')
        plt.plot(readout_point,vis, 'ro')
        # print(f'VIS={vis*100}%, vP1= {readout_point} mV')
    return readout_point, vis