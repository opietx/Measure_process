# Q1
# uuids = [1655734070027283874,1655734380021283874,1655734660764283874,1655734380021283874,1655734743606283874,1655735008004283874]
# title = 'base, reverse Rabi'

uuids = [1655735686203283874,1655735962916283874,1655735962916283874,1655736118072283874,1655736270386283874,1655736571903283874]
title = 'base, forward Rabi'

# uuids = [1655738282780283874,1655738366918283874,1655738408243283874,1655738493050283874,1655738628892283874,1655738898616283874]
# title = '200mK, reverse Rabi'

# uuids = [1655738282780283874,1655738366918283874,1655738408243283874,1655738493050283874,1655738628892283874,1655738898616283874]
# title = '200mK, forward Rabi'



qubit = 1
dic_data = {'Rabi time':None}
times = [0,100,500,1000,2000,4000]

for ii, uuid  in enumerate(uuids):
    ds = load_by_uuid(uuid)
    y = ds['read1'].y()
    dic_data[str(times[ii])+' us'] = y
    

    

dic_data['Rabi time'] = ds['read1'].x()
plt.figure()
for ii,key in enumerate(times):
    plt.plot(dic_data['Rabi time'],dic_data[str(times[ii])+' us'],'o-',label=key)
plt.ylabel('Spin fraction')
plt.xlabel('Rabi_time [ns]')
plt.legend()
plt.title(title)
plt.show()     


# data_dump(dic_data, 'Rabi_wait_time','12mK_forward')


#%%
from measure_process.fits.Rabi import  fit_Rabi, Rabi_formula

ds = load_by_uuid(1655822946129283874)
y = ds['read1'].y()
x = ds['read1'].x()

a=fit_Rabi(x*1e-9,y,'aa',True)

#%%
uuids = [1655797073603283874,1655798006137283874, 1655798565278283874]
lables=['rough-calib','ramsey-calib','rough-calib2']
plt.figure()
for ii, uuid  in enumerate(uuids):
    ds = load_by_uuid(uuid)
    x = ds['read1'].x()
    y = ds['read1'].y()
    plt.plot(x,y,'o-', label=lables[ii])
plt.ylabel('Spin fraction')
plt.xlabel('Rabi_time [ns]')
plt.legend()
plt.show()