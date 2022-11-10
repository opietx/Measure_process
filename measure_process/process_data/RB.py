import matplotlib.pyplot as plt
import numpy as np

from measure_process.Init_file import data_dump
from core_tools.data.ds.data_set import load_by_uuid
from measure_process.fits.RB import fit_RB



qubit = 2
plot = True 
uuid = 1668010887753974087
ds = load_by_uuid(uuid)
# view = data_plotter(ds)

ds_avg = ds[f'read{qubit}'][:].average('x')
params, errors = fit_RB(ds_avg.x()[:], ds_avg.y()[:], plot=plot)
