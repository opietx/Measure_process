import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid
from measure_process.fits.RB import fit_RB

#%%

qubit = 1   
plot = True 
uuid = 1668192989006974076
ds = load_by_uuid(uuid)
# view = data_plotter(ds)
ds_avg = ds[f'read{qubit}'][:].average('x')
x = (ds_avg.x()[:]).astype(int)-1# add this for old measurements where we had 1.1 number of cliffords like-array
params1, errors = fit_RB(x, ds_avg.y()[:], plot=plot, name='Q2_X')
