import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid
from measure_process.fits.T1 import  Decay_formula, fit_decay_T1

#%% 1 dataset averaged all
qubit = 2
plot = True 
uuid = 1658793092358974076      
ds = load_by_uuid(uuid)
# data_plotter(ds)
ds_avg = ds[f'read{qubit}'].average('x')
x,y = 1e-9*ds_avg.x(), ds_avg.y()
#it is 1-y for T1 and only y for T1_PSB
params, errors = fit_decay_T1(x,y,uuid, plot=plot)
print(errors['T1'])

#%%

# qubit = 1
# plot = True 
# uuid = 1659464931048974087     
# ds = load_by_uuid(uuid)
# data_plotter(ds)
# ds_avg = ds[f'read{qubit}'].average('x')

# temp=200
# # name = f'Q{qubit}'
# name = f'Q{qubit}_baseline'

# dic_data = {'time':x,f'Q{qubit}':y}
# data_dump(dic_data,f'T1/{temp}mK',name)



#%% to select multiple datasets parts
y_multiple = []

qubit = 1
plot = True 
uuid = 1659464931048974087    
ds = load_by_uuid(uuid)
# data_plotter(ds)


for ii in range(1,3):
    y_multiple.append(ds[f'read{qubit}'][ii].y())

y_avg = sum(y_multiple)/len(y_multiple)


#it is 1-y for T1 and only y for T1_PSB
params, errors = fit_decay_T1(1e-9*ds[f'read{qubit}'][0].x(), y_avg, uuid, plot=plot)
print(errors['T1'])


#%% 1Q Load data, multiple datasets
temps = [200,250,300,350,400,450,500,550,600,650,700,750,800,850]



# #Q1
# qubit = 1
# uuids = [1658446151857974076 ,1658461592095974076 ,1658477032456974076 ,1658492471876974076 ,1658507899401974076 ,1658523333405974076 ,1658538768833974076 ,
#          1658566513665974076 ,1658573671096974076 ,1658632283523974076 ,1658645223592974076 ,1658724765041974076 ,1658709404104974076 ,1658681273214974076]


#Q2
# qubit = 2
# uuids = [1658448580282974076 ,1658464020984974076 ,1658479461667974076 ,1658494900575974076 ,1658510327909974076 ,1658525761833974076 ,1658541197193974076 ,
#          1658567544112974076 ,1658574692842974076  ,1658634380914974076 ,1658647319994974076 ,1658727232464974076 ,1658711872298974076 ,1658683101392974076]

# qubit = 5
# uuids = [1658451007116974076 ,1658466448065974076 ,1658481887796974076 ,1658497327738974076 ,1658512754015974076 ,1658528188199974076 ,1658543623778974076 ,
#          1658568563465974076 ,1658575712217974076 ,1658636474971974076 ,1658649414132974076 ,1658662336302974076 ,1658714337621974076 ,1658698993512974076]


qubit = 6
uuids = [1658453436606974076 ,1658468877129974076 ,1658484316924974076 ,1658499757821974076 ,1658515182417974076 ,1658530616363974076 ,1658546052137974076 ,
         1658569593898974076 ,1658576732890974076 ,1658638568750974076 ,1658651508171974076 ,1658732159113974076 ,1658716798852974076 ,1658701456105974076]

#%%##
dic_data = {'Temp':temps,f'Q{qubit}':np.empty(len(temps)),'low':np.empty(len(temps)),'top':np.empty(len(temps))}
for ii, uuid  in enumerate(uuids):
    ds = load_by_uuid(uuid)
    ds_avg = ds[f'read{qubit}'].average('x')
    params,errors = fit_decay_T1(1e-9*ds_avg.x(), ds_avg.y(), uuid, plot=False)
    dic_data[f'Q{qubit}'][ii] = params['T1'].value
    dic_data['low'][ii] = errors['T1'][0]
    dic_data['top'][ii] = errors['T1'][1]
    


plt.figure()
plt.errorbar(np.array(temps),dic_data[f'Q{qubit}'], yerr=dic_data['low'],fmt='o-',label=qubit)
plt.xlabel('T_mc [mK]')
plt.ylabel('T1 [s]')
plt.legend()
plt.show()        


#dump data
data_dump(dic_data,'T1',f'Q{qubit}')