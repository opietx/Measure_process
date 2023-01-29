import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid
from measure_process.fits.T1 import  Decay_formula, fit_decay_T1
#%%

qubit = 5
plot = True 
uuid = 1658795515753974076   
ds = load_by_uuid(uuid)
# data_plotter(ds)
ds_avg = ds[f'read{qubit}'].average('x')
# ds_avg_y = ds_avg.y()-ds_avg.y().min()
# ds_avg_y /= ds_avg_y.max()
#it is 1-y for T1 and only y for T1_PSB
params, errors = fit_mixed_decay_T1(1e-9*ds_avg.x(), 1-ds_avg.y(), uuid, plot=plot)
