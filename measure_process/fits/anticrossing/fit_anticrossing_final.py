import time
import numpy as np
import matplotlib.pyplot as pt

from core_tools.data.SQL.queries.dataset_gui_queries import query_for_measurement_results
from measure_process.fits.anticrossing.ds2array import ds2xarray

from projects.xarray_utils.data import xy_coord_values
from projects.xarray_utils.vector_fields import create_vf, vf2polar, vf_hsv_plot, VectorField

from measure_process.fits.anticrossing.images.filters import gradient_vf, subtract_noise, compress
from measure_process.fits.anticrossing.images.anticrossing import plot_cross, fit_anti_crossing, estimate_center
from core_tools.data.ds.data_set import load_by_uuid

def plot_add_text(msgs=None, warning=None):
    ax = pt.gca()
    if msgs:
        y = 0
        for msg in msgs:
            pt.text(0.5, y, msg,
                    transform = ax.transAxes,
                    color='cyan',
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    fontsize=14)
            y += 0.06
    if warning:
            pt.text(0.5, 0.98, f'*** {warning} ***',
                    transform = ax.transAxes,
                    color='red',
                    horizontalalignment='center',
                    verticalalignment='top',
                    fontsize=16)


def plot_data(data, name, title, points=None, msgs=None, warning=None, **kwargs):
    pt.figure()
    data.plot(**kwargs)
    if points is not None:
        plot_cross(points)
    plot_add_text(msgs, warning)
    pt.title(title)
    # pt.savefig(f'{directory}/{name}.png')
    # pt.close()
    pt.show()


def plot_hsv(xrdata, name, title, points=None, msgs=None, warning=None, ):
    pt.figure(figsize=(5.6, 4.8))
    vf_hsv_plot(xrdata)
    if points is not None:
        plot_cross(points)
    plot_add_text(msgs, warning)
    pt.title(title)
    # pt.savefig(f'{directory}/{name}.png')
    # pt.close()
    pt.show()



def fit_anticrossing(uuid,config, averaged=True, plot=True,verbose=False):
        
    directory = 'figures_ac'
    
    
    # parameters for the anti-crossing
    vf,sigma_gradient,sigma_estimate,noise_subtraction,estimate_start,all_legs,leg_length,sigma_model,angle_center = config
    legs = ['c',1,2,3,4] if all_legs else  ['c',1,3]

    
    
    # disable automatic display of plots
    pt.ioff()
    
    d = load_by_uuid(uuid)
    
    
    title = f'{d.run_timestamp:%Y-%m-%d %H.%M.%S} ({d.exp_uuid})'
    pt_title = title
    dt = f'{d.run_timestamp:%Y-%m-%d %H.%M.%S}'
    
    dsx = ds2xarray(d)
    
    # skip 1D scans
    
    # all scans use channel 4 of digitizer
    ch4 = dsx['ch3']
    if plot:
        plot_data(ch4, f'ch4_{dt}', pt_title)
    
    ch4_grad = gradient_vf(ch4, sigma=sigma_estimate)
    ch4_mag, ch4_angle = vf2polar(ch4_grad)
    
    start_x, start_y, ac_img = estimate_center(ch4_mag, ch4_angle)
    
    if 1:
        if plot:
            plot_data(ac_img, f'ch4_{dt}_anticrossing', pt_title)
    
    ch4_grad = gradient_vf(ch4, sigma=sigma_gradient)
    ch4_mag, ch4_angle = vf2polar(ch4_grad)
    if plot:
        plot_data(ch4_mag, f'ch4_{dt}_mag', pt_title)
    
    if noise_subtraction:
        # about 15% of image is honeycomb gradient. Rest is 'noise'
        # NOTE: subtracting noise doesn't seem to improve the fit quality.
        ch4_mag = subtract_noise(ch4_mag, percentile=85)
    
    #    plot(ch4_grad.sel(XY=0), f'ch4_{dt}_x', pt_title)
    #    plot(ch4_grad.sel(XY=1), f'ch4_{dt}_y', pt_title)
    #    plot(ch4_phase, f'ch4_{nr}_phase', pt_title, cmap='hsv')
    
    if not estimate_start:
        start_x, start_y = 1693,1563
    
    if all_legs:
        fit_param_names = ['center_x', 'center_y', 'width', 'angle_right', 'angle_left']
    else:
        fit_param_names = ['center_x', 'center_y', 'width', 'angle_left']
    
    # do the actual fitting
    start = time.perf_counter()
    img = ch4_grad if vf else ch4_mag
    # print( start_x, start_y)
    fitted_cross = fit_anti_crossing(img, fit_param_names, sigma_model,
                                     legs, leg_length, start_x, start_y, vf,
                                     fit_xy_rotation=0.01,
                                     center_angle=angle_center)
    duration = (time.perf_counter() - start)
    
    ### extract data from result
    attrs = fitted_cross.attrs
    match = attrs['match']
    match_legs = attrs['match_legs']
    fp = attrs['fitted_model']
    center_x, center_y = attrs['center']
    width = attrs['width']
    angle_right = fp['angle_right'] if all_legs else  0.0
    angle_left = fp['angle_left']
    points = attrs['points']
    
    x,y = xy_coord_values(fitted_cross)
    
    ### determine whether fit is of sufficient quality
    ok = (
        attrs['success']
        and (
             (attrs['match'] > 0.8 and width > 1.5)
             or (match_legs[1] > 0.2 and match_legs[3] > 0.2 and width > 2.0)
             or (match_legs[2] > 0.2 and match_legs[4] > 0.2 and width > 2.0)
             or (match_legs[1] > 0.15 and match_legs[3] > 0.15
                 and match_legs[2] > 0.15 and match_legs[4] > 0.15)
             )
        )
    if not ok:
        print(f'UUID:{uuid}  *** FAILED' )
    
    
    ### print and plot the results
    warning = None
    if verbose:
        print(title)
        print(f'  duration {duration:5.2f} s')
        print(f'  start:  {start_x:.1f}, {start_y:.1f}')
        print(f'  center: {center_x:.1f},  {center_y:.1f}; width: {width:.1f} mV')
        print(f'  angles: {angle_left*180/np.pi:.1f}, {angle_right*180/np.pi:.1f} degrees')
        print(f'  Match {match:3.0%}  ({match_legs[1]*100:2.0f} {match_legs[2]*100:2.0f} '
              f'{match_legs[3]*100:2.0f} {match_legs[4]*100:2.0f})  {width:3.1f} mV')
    
        
        if not ok:
            warning = 'FAILED'
            print(f'  *** FAILED *** ({attrs["success"]})' )
        elif abs(center_y)-width/2 > 0:
            warning = 'NOT CENTERED'
            print('  *** NOT CENTERED ***' )
    
    if plot:
        pt_title = title + f' ({match:5.1%})'
        
        msgs = [
            f'{center_x:3.1f},{center_y:3.1f} mV ({width:3.1f} mV)',
            f'[{match_legs[1]*100:2.0f} {match_legs[2]*100:2.0f} {match_legs[3]*100:2.0f} {match_legs[4]*100:2.0f}]'
            ]
        
        plot_data(ch4, f'ch4_{dt}_fit.png', pt_title, points, msgs, warning)
        
        if not vf:
            if 1:
                ch4_mag = compress(ch4_mag)
            plot_data(ch4_mag, f'ch4_{dt}_mag', pt_title)
            plot_data(ch4_mag, f'ch4_{dt}_mag_f', pt_title, points, msgs, warning)
        
        if vf:
            if 1:
                ch4_mag = compress(ch4_mag, subtract_low=True)
            ch4_rphi2 = create_vf([ch4_mag, ch4_angle], VectorField.RPhi)
            plot_hsv(ch4_rphi2, f'ch4_{dt}_grad', pt_title)
            plot_hsv(ch4_rphi2, f'ch4_{dt}_grad_f', pt_title, points, msgs, warning)
        
        if 0:
            # plot generated anti-crossing
            if vf:
                plot_hsv(fitted_cross, f'ch4_{dt}_cross_vf', pt_title)
            else:
                plot_data(fitted_cross, f'ch4_{dt}_cross', pt_title)
    plot_data(ch4_mag, f'ch4_{dt}_mag_f', pt_title, points, '', None)
    return points
        
