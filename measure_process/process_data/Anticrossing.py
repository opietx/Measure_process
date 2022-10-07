import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from measure_process.fits.anticrossing.fit_anticrossing_final import  fit_anticrossing

#%%

config_params = [False,2.5,0.1,True,True,True,0.5,7,np.deg2rad(55)]


ch='ch3'
plot = True 
uuid = 1662057789052974076    
params= fit_anticrossing(uuid,config_params,averaged=True, plot=plot,verbose=False)
point1 = params[0]

#%%
temps = [25,50,75,100,125,150,175,200,225,250]



uuids = [1662030571359974076 ,1662033788267974076 ,1662037239069974076 ,1662040682483974076 ,1662044176659974076 ,1662047639350974076 ,1662051104532974076 
         ,1662054427599974076 ,1662057789052974076 ,1662061101022974076 ]

#%%

config_params = [False,2.5,0.1,True,True,True,0.5,7,np.deg2rad(55)]
dic_data = {'Temp':temps,'x1':np.empty(len(temps)),'x2':np.empty(len(temps)),'y1':np.empty(len(temps)),'y2':np.empty(len(temps))}
points = []
for ii, uuid  in enumerate(uuids):
    print(f'temp step: {temps[ii]}mK')
    params= fit_anticrossing(uuid,config_params,averaged=True, plot=False,verbose=False)
    
    dic_data['x1'][ii] = params[0][0]
    dic_data['x2'][ii] = params[1][0]
    dic_data['y1'][ii] = params[0][1]
    dic_data['y2'][ii] = params[1][1]
    

    




plt.figure()
# plt.errorbar(np.array(temps),dic_data[f'{sensor}'], yerr=dic_data['low'],fmt='o-',label=sensor)
plt.plot(temps,dic_data['x1']/dic_data['x1'][0],'-o')
plt.plot(temps,dic_data['y1']/dic_data['y1'][0],'-o')
plt.plot(temps,dic_data['x2']/dic_data['x2'][0],'-o')
plt.plot(temps,dic_data['y2']/dic_data['y2'][0],'-o')

plt.xlabel('T_mc [mK]')
plt.ylabel('position')
# plt.legend()
plt.show()        

data_dump(dic_data,'Stability','test')