import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid
from measure_process.fits.T1 import  Decay_formula, fit_decay_T1

#%% 1 dataset averaged all
qubit = 56
plot = True 
uuid = 1661259476593974087      
ds = load_by_uuid(uuid)
data_plotter(ds)
ds_avg = ds[f'read{qubit}_no_ramp'].average('x')

x,y = 1e-9*ds_avg.x(), ds_avg.y()

params, errors = fit_decay_T1(x,y, uuid, plot=plot)
print(errors['T1'])

#%% Dump 1dataset data
temp=800
state = 'init'
dic_data = {'time':x,f'Q{qubit}':y}
data_dump(dic_data,f'T1_PSB/SD{qubit}/{temp}mK',f'{state}')


#%% to select multiple datasets parts
y_multiple = []

qubit = 56
plot = True 
uuid = 1661259049716974087    
ds = load_by_uuid(uuid)
# data_plotter(ds)


for ii in range(1,3):
    y_multiple.append(ds[f'read{qubit}_no_ramp'][ii].y())

y_avg = sum(y_multiple)/len(y_multiple)
x,y = 1e-9*ds[f'read{qubit}_no_ramp'][0].x(), y_avg

#it is 1-y for T1 and only y for T1_PSB
params, errors = fit_decay_T1(x,y, uuid, plot=plot)
print(errors['T1'])


#%% 1Q Load data, multiple datasets
temps = [200,250,300,350,400,450,500,550,600,650,700,750,800,850,900]
# SD=12
#first set
# state = 'T_dd'
# uuids = [1658951492992974076 ,1658953217163974076 ,1658954729351974076 ,1658956361528974076 ,1658957995917974076 ,1658959586500974076 ,1658960881591974076 ,
#           1658962177406974076 ,1658963472453974076 ,1658964741305974076 ,1658965918095974076 ,1658967056334974076 ,1658968078517974076 ,1658969103578974076 ]


#first set
# state = 'T_dd'
# uuids = [1659291500175974087 ,1659293044041974087 ,1659294340591974087 ,1659295759434974087 ,1659297194299974087 ,1659298445415974087 ,1659299368800974087 ,
#          1659300281387974087 ,1659301195316974087 ,1659302063591974087 ,1659302794274974087 ,1659303547290974087 ]

# state= 'T_uu'
# uuids= [1659291878884974087 ,1659293414927974087 ,1659294712663974087 ,1659296131350974087 ,1659297563390974087 ,1659298654315974087 ,1659299581465974087 ,
#         1659300493723974087 ,1659301408375974087 ,1659302206876974087 ,1659302942083974087 ,1659303694980974087]




#after retunning
# SD=56
# state='T_dd'
# uuids = [1661258436010974087 ,1661259049716974087 ,1661259476593974087 ,1661260021305974087 ,1661260568525974087 ,1661261165923974087 ,1661261667234974087 ,
#           1661262145603974087 ,1661262608660974087 ,1661263054881974087 ,1661263479504974087 ,1661263935944974087 ,1661264367598974087 ,1661264753583974087 ,
#           1661265142085974087 ]

# state = 'T_uu'
# uuids = [1661258531203974087 ,1661259144639974087 ,1661259571630974087 ,1661260115225974087 ,1661260663119974087 ,1661261247991974087 ,1661261750051974087 ,
#           1661262228254974087 ,1661262673090974087 ,1661263119198974087 ,1661263544434974087 ,1661264000458974087 ,1661264420392974087 ,1661264806147974087 ,
#           1661265194756974087 ]

#%%##
dic_data = {'Temp':temps,f'{state}':np.empty(len(temps)),'low':np.empty(len(temps)),'top':np.empty(len(temps))}
for ii, uuid  in enumerate(uuids):
    ds = load_by_uuid(uuid)
    ds_avg = ds[f'read{SD}_no_ramp'].average('x')
    params,errors = fit_decay_T1(1e-9*ds_avg.x()[:], ds_avg.y()[:], uuid, plot=False)
    dic_data[f'{state}'][ii] = params['T1'].value
    dic_data['low'][ii] = errors['T1'][0]
    dic_data['top'][ii] = errors['T1'][1]
    


plt.figure()
plt.errorbar(np.array(temps),dic_data[f'{state}'], yerr=dic_data['low'],fmt='o-',label=f'SD:{SD} {state}')
plt.xlabel('T_mc [mK]')
plt.ylabel('T1 [s]')
plt.legend()
plt.show()        


#dump data
data_dump(dic_data,f'T1_PSB/SD{SD}_before_tune',f'{state}')