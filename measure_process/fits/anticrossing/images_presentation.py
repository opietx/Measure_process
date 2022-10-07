import time
import numpy as np
import matplotlib.pyplot as pt

from core_tools.data.ds.data_set import load_by_uuid
from core_tools.data.ds.ds2xarray import ds2xarray
from core_tools.data.SQL.queries.dataset_gui_queries import query_for_measurement_results

from projects.xarray_utils.vector_fields import create_vf, vf2polar, vf_hsv_plot, VectorField

from images.filters import gradient_vf, subtract_noise, gradient_sobel, compress

from projects.xarray_utils.data import derive

import sds_test_database
#sds_test_database.setup_local_and_remote_other('veldhorst_data')
sds_test_database.setup_local_other('veldhorst_data')



def plot(data, name, title, **kwargs):
    pt.figure()
    data.plot(**kwargs)
    pt.title(title)
    pt.savefig(f'{directory}/{name}.png')
    pt.close()


def plot_hsv(xrdata, name, title):
    pt.figure(figsize=(5.6, 4.8))
    vf_hsv_plot(xrdata)
    pt.title(title)
    pt.savefig(f'{directory}/{name}.png')
    pt.close()


res = query_for_measurement_results.search_query(
        start_time='2021-09-06',
        end_time='2021-09-08 23:59',
        name='0D_digitizer_',
#        remote=True
        )


directory = 'figures_plots'


# parameters for the anti-crossing
vf = True
sigma_gradient = 1.3
noise_subtraction = False

# disable automatic display of plots
pt.ioff()

for e in res:
    uuid = e.uuid

    d = load_by_uuid(uuid, copy2localdb=True)

    title = f'{d.run_timestamp:%Y-%m-%d %H.%M.%S} ({d.exp_uuid})'
    pt_title = title
    dt = f'{d.run_timestamp:%Y-%m-%d %H.%M.%S}'

    dsx = ds2xarray(d)

    # skip 1D scans
    if len(dsx.dims) != 2:
        continue

    print(pt_title)
    # all scans use channel 4 of digitizer
    ch4 = dsx['ch4']
    plot(ch4, f'ch4_{dt}__', pt_title)

    ch4__c = derive(ch4)
    ch4__c.data = compress(ch4, upper=99.5, lower=0.5)
    plot(ch4__c, f'ch4_{dt}__x', pt_title)

    sobel = gradient_sobel(ch4)
    sobel_mag, sobel_angle = vf2polar(sobel)
    plot(sobel_mag, f'ch4_{dt}_sobel_mag', pt_title)
    sobel_rphi = create_vf([sobel_mag, sobel_angle], VectorField.RPhi)
    plot_hsv(sobel_rphi, f'ch4_{dt}_sobel', pt_title)

    sobel_mag_x = compress(sobel_mag)
    plot(sobel_mag_x, f'ch4_{dt}_sobel_mag_x', pt_title)
    sobel_mag_x = compress(sobel_mag, subtract_low=True)
    sobel_rphi = create_vf([sobel_mag_x, sobel_angle], VectorField.RPhi)
    plot_hsv(sobel_rphi, f'ch4_{dt}_sobel__x', pt_title)

#    ch4_dx = derive(ch4)
#    ch4_dx.data[:,1:] = ch4.values[:,:-1] - ch4.values[:,1:]
#    ch4_dx.data[:,0] = 0
#    ch4_dy = derive(ch4)
#    ch4_dy.data[1:,:] = ch4.values[:-1,:] - ch4.values[1:,:]
#    ch4_dy.data[0,:] = 0
#    ch4_dif_mag = derive(ch4)
#    ch4_dif_mag.data = np.sqrt(ch4_dx.values**2 + ch4_dy.values**2)
#    plot(ch4_dif_mag, f'ch4_{dt}__dif_mag', pt_title)

    ch4_grad = gradient_vf(ch4, sigma=0.9)
    ch4_mag, ch4_angle = vf2polar(ch4_grad)
    ch4_rphi2 = create_vf([ch4_mag, ch4_angle], VectorField.RPhi)
    plot_hsv(ch4_rphi2, f'ch4_{dt}_grad__', pt_title)
    plot(ch4_mag, f'ch4_{dt}_mag__', pt_title)

    ch4_mag_x = compress(ch4_mag)
    plot(ch4_mag_x, f'ch4_{dt}_mag__x', pt_title)
    ch4_mag_x = compress(ch4_mag, subtract_low=True)
    ch4_rphi2_x = create_vf([ch4_mag_x, ch4_angle], VectorField.RPhi)
    plot_hsv(ch4_rphi2_x, f'ch4_{dt}_grad__x', pt_title)

    ch4_grad = gradient_vf(ch4, sigma=sigma_gradient)
    ch4_mag, ch4_angle = vf2polar(ch4_grad)
    plot(ch4_mag, f'ch4_{dt}_z_mag', pt_title)

    ch4_rphi2 = create_vf([ch4_mag, ch4_angle], VectorField.RPhi)
    plot_hsv(ch4_rphi2, f'ch4_{dt}_z_grad', pt_title)

    if noise_subtraction:
        # about 15% of image is honeycomb gradient. Rest is 'noise'
        ch4_mag = subtract_noise(ch4_mag, percentile=85)

#    plot(ch4_grad.sel(XY=0), f'ch4_{dt}_x', pt_title)
#    plot(ch4_grad.sel(XY=1), f'ch4_{dt}_y', pt_title)
#    plot(ch4_phase, f'ch4_{nr}_phase', pt_title, cmap='hsv')




