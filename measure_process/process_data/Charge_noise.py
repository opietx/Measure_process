import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid
from measure_process.fits.Charge_noise import  fit_PSD

import numpy as np
import scipy as sp
from core_tools.data.ds.data_set import load_by_uuid


def PSD(x, fs, nperseg):
    n_iter = int(x.size/nperseg)
    
    freq = np.fft.rfftfreq(nperseg, d=1/fs)
    FT_coeff = np.zeros(freq.shape)
    for i in range(n_iter):
        FT_coeff += 2*np.abs(np.fft.rfft(x[i*nperseg:(i+1)*nperseg]))**2/fs/nperseg
    
    return freq, FT_coeff/n_iter

def to_PSD(uuid, n_seg = 2, lever_arm=0.18, new_slope =1, old_slope=1):
    '''
    if old_slope is given is to be removed from measurement to be replaced with a more accurate version
    '''
    ds_1 = load_by_uuid(uuid)
    y = ds_1.m1()*old_slope/new_slope
    
    
    
    
    f_s= 1/ds_1.m1.x()[1]
    nperseg = 2**int(np.log2(int(y.size/n_seg)))
    
    freq, Sxx_SD = PSD(y, fs =f_s, nperseg = nperseg)#makes the PSD of the signal
    Sxx = Sxx_SD #the background is not removed anymore
    Sxx *= lever_arm*lever_arm
    
    return freq, Sxx

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def obtain_slope_give_voltage(uuid,meas_point):
    from scipy.signal import savgol_filter
    """ Allows the user to click three points. This method will calculate the (positive)
        slope between point 1 and 2 and defines the coordinates of point 3.
        Args:
            ds (dataset) : dataset with charge stability diagram and gate voltage in mV
            background_subtraction (bool) : include background subtraction as clickable point
        Returns:
            slope (double) : the slope between point 1 and point 2 (always positive)
            noise_meas (double) : the x coordinate of the noise measurement
            bg_meas (double) : the x coordinate of the background measurement
    """
    ds = load_by_uuid(uuid)
    x = ds.m1.x()[20:-20]
    y = ds.m1.y()[20:-20]
    
    
    y_dif = np.gradient(y)
    x_dif = x
    y_dif = savgol_filter(y_dif, 51, 3) # window size 51, polynomial order 3
    y_dif_search = y_dif[:int(len(y_dif)/2)]
    
    arg_max = find_nearest(x,meas_point)
    
    slope = y_dif[arg_max]*1000 #gate in mV
    
    
    
    return slope

#%%
plot = True 
uuid =1659360940238974087#     
new_slope =obtain_slope_give_voltage(1659360819650974087, 1534.9369138273353)
x, y = to_PSD(uuid, n_seg=20, lever_arm=0.185, new_slope=new_slope, old_slope=6.602043148313142e-12)    

params, errors, [x,y] = fit_PSD(x,y,name=uuid,plot=plot)
x,y = x[:500], y[:500]
dic_psd = {'Freq':x,f'PSD':y}
data_dump(dic_psd,'Charge_noise',f'500mk') #if you want to export the PSD aswell

#%%
temps = [200,300,400,500,600,700,800,900,950,1000,1050,1100,1200,1300]

# sensor = 'SD1'
# lever_arm = 0.185
# uuids = [1659354850881974087 ,1659357709833974087 ,1659359708963974087 ,1659360940238974087 ,1659362851218974087 ,1659363786921974087 ,1659366053628974087 ,
#           1659368874303974087 ,1659427856696974087 ,1659371323853974087 ,1659429945815974087 ,1659430881193974087 ,1659433832236974087 ,1659439070335974087 ]

# slopes = [1.4299499961589247e-11,1.167786788791675e-11 ,8.225442572684299e-12 ,6.602043148313142e-12 ,5.43122026258292e-12 ,4.356922431141352e-12 ,3.6302939682948103e-12 ,
#           3.6566233537955345e-12 ,2.9800366147027407e-12 ,2.9949990397863812e-12 ,2.6185045321909965e-12 ,2.382553761033205e-12 ,2.335449840442744e-12 ,2.164297172142734e-12 ]

# meas_points = [1534.995030060120 ,1534.9950300599844 ,1534.8206813624058 ,1534.9369138273353 ,1534.9659719433362 ,1534.5882164323136 ,1534.7044488971046 ,1534.5010420834762 ,
#                 1534.4971869683432 ,1534.2976352696023 ,1534.758710014321 ,1534.49718696822 ,1534.2066058057 ,1534.613419433029 ]

# uuids_CP = [1659354729808974087 ,1659357589274974087 ,1659359589008974087 ,1659360819650974087 ,1659362731340974087 ,1659363667032974087 ,1659365933338974087 ,1659368754515974087 ,
#             1659427736597974087 ,1659371204071974087 ,1659429825838974087 ,1659430761400974087 ,1659433711706974087 ,1659438949880974087 ]


sensor = 'SD2'
lever_arm = 0.19
uuids = [1659356478938974087 ,1659358312307974087 ,1659359150897974087 ,1659361476491974087 ,1659362250391974087 ,1659364406830974087 ,1659365456716974087 ,
          1659369410168974087 ,1659428375722974087 ,1659370816981974087 ,1659429429006974087 ,1659431957061974087 ,1659433316827974087 ,1659439582651974087 ]

slopes = [1.3562064590959033e-11 ,1.175648791337511e-11 ,9.972096278564048e-12 ,9.212576078496227e-12 ,7.610418803063019e-12 ,6.365134997642369e-12 ,6.1328745276286226e-12,
          6.2255371146159394e-12 ,2.9800366147027407e-12 ,5.380689665014759e-12 ,3.846599822401599e-12 ,4.030446492001069e-12 ,3.735851631626064e-12 ,3.3543491133145807e-12 ]

meas_points = [1524.351988459782 ,1524.3519884597329 ,1524.2357559948023 ,1524.1776397623175 ,1524.1485816460845 ,1523.9161167161853 ,1523.654593670093 ,1523.6836517862334 ,
                1534.4971869683432 ,1523.8580004836278 ,1523.480244972932 ,1523.65459367044 ,1523.6545936704456 ,1523.65459367 ]

uuids_CP = [1659356358011974087 ,1659358192096974087 ,1659359030833974087 ,1659361356730974087 ,1659362130263974087 ,1659364286970974087 ,1659365336661974087 ,1659369290162974087 ,
            1659428255305974087 ,1659370696310974087 ,1659429308382974087 ,1659431836856974087 ,1659433197063974087 ,1659439462506974087 ]
#%%
dic_data = {'Temp':temps,f'{sensor}':np.empty(len(temps)),'low':np.empty(len(temps)),'top':np.empty(len(temps))}
for ii, uuid  in enumerate(uuids):
    new_slope = obtain_slope_give_voltage(uuids_CP[ii], meas_points[ii])
    x, y = to_PSD(uuid, n_seg=20, lever_arm=lever_arm, new_slope=new_slope, old_slope=slopes[ii])    
    params, errors = fit_PSD(x,y,name=uuid,plot=True)
    # dic_psd = {'Freq':x,f'{temps[ii]}mK':y}
    # data_dump(dic_psd,'Charge_noise',f'{sensor}_{temps[ii]}mK') #if you want to export the PSD aswell
    
    # extract amplitude 
    # name = f'{sensor}'
    # dic_data[f'{sensor}'][ii] = np.sqrt(abs(params['P_noise'].value))
    # dic_data['low'][ii] = 0.5*(errors['P_noise'][0]/abs(params['P_noise'].value))*np.sqrt(abs(params['P_noise'].value))
    # dic_data['top'][ii] = 0.5*(errors['P_noise'][1]/abs(params['P_noise'].value))*np.sqrt(abs(params['P_noise'].value))
    
    #extract alpha
    name = f'alpha_{sensor}'
    dic_data[f'{sensor}'][ii] = abs(params['alpha'].value)
    dic_data['low'][ii] = errors['alpha'][0]
    dic_data['top'][ii] = errors['alpha'][1]  


plt.figure()
plt.errorbar(np.array(temps),dic_data[f'{sensor}'], yerr=dic_data['low'],fmt='o-',label=sensor)
plt.xlabel('T_mc [mK]')
plt.ylabel('P_noise')
plt.legend()
plt.show()        

data_dump(dic_data,'Charge_noise',name)

#dump data

#%% extract coulomb peaks to a file for origin
for ii, uuid in enumerate(uuids_CP):
    name = str(temps[ii])
    ds = load_by_uuid(uuid)
    x,y = ds.m1.x(), ds.m1.y()

    
    dic_data = {f'{sensor}':x,f'I, {name}mK':y}
    data_dump(dic_data,f'Charge_noise/{sensor}',f'{name}mK')


