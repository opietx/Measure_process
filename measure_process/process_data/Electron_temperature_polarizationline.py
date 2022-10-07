import numpy as np
import matplotlib.pyplot as plt

from measure_process.fits.Polarization_lines import  Fermi_center
from core_tools.data.ds.data_set import load_by_uuid

#%%

ds = load_by_uuid(1663155505768974076)

ds_avg = ds['ch3'].average('x')
y =  ds_avg.y()-(ds_avg.y()).min()
y = y/y.max()
a= Fermi_center(ds_avg.x(),y, plot=True)

#%%
uuids = [1663152483615974076,1663152820344974076,1663152946284974076,1663153042100974076,1663153438991974076,1663153835903974076,
         1663154172551974076,1663154388859974076,1663154454545974076,1663154610648974076,1663154796815974076,1663155073271974076,
         1663155349687974076,1663155505768974076]
temps = [25,50,75,100,125,150,175,200,225,250,300,350,400,500]

#%%
temps = np.arange(30,210,10)
temps = np.append(temps,np.array([400,500]))
uuids = [1663159667259974076,1663159943750974076,1663160250338974076,1663160376291974076,1663160472005974076,1663160567876974076,1663160964822974076,1663161361741974076,
         1663161487713974076,1663161583525974076,1663161950212974076,1663162076088974076,1663162171742974076,1663162508348974076,1663162634226974076,1663162730078974076,
         1663163066845974076,1663163253094974076,
         1663155349687974076,1663155505768974076]#400,500

#%%
dic_data = {'Temp':np.array(temps)/1000,f'Te':np.empty(len(temps))}

for ii, uuid  in enumerate(uuids):
    ds = load_by_uuid(uuid)
    
    ds_avg = ds['ch3'].average('x')
    y =  ds_avg.y()-(ds_avg.y()).min()
    y = y/y.max()
    a= Fermi_center(ds_avg.x(),y, plot=True)
    dic_data[f'Te'][ii] = a[2]
    

plt.figure()
plt.plot(np.array(temps),dic_data[f'Te'])
plt.xlabel('T_mc [mK]')
plt.ylabel('Te [mK]')
plt.show()        


#dump data
# data_dump(dic_data,'Electron_temp',f'Te')

    