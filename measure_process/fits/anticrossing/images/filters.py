from numbers import Number
import numpy as np
import scipy.ndimage as ndimage

from projects.xarray_utils.vector_fields import create_vf, VectorField
from projects.xarray_utils.data import derive

def gradient_magnitude(data_array, sigma=1.0):
    res = derive(
            data_array,
            f'gradient(sigma={sigma})',
            )
    res.data = ndimage.gaussian_gradient_magnitude(res.data, sigma=sigma, mode='nearest')
    return res

def gradient_vf(data, sigma=1.0):
    '''
    First dimension is 'y'. Second dimension is 'x'.

    Returns: 3D DataArray with dim 'yx'; 'yx'==0 is d/dy, 'yx'==1 is d/dx.
    '''
    x = derive(data, f'gradient(sigma={sigma})')
    y = derive(data, f'gradient(sigma={sigma})')

    x.data = ndimage.gaussian_filter(data, order=(0,1), sigma=sigma, mode='nearest')
    y.data = ndimage.gaussian_filter(data, order=(1,0), sigma=sigma, mode='nearest')

    return create_vf([x,y], VectorField.XY)

def gradient_sobel(data):
    sobel_x = derive(data, 'sobel()')
    sobel_y = derive(data, 'sobel()')
    sobel_x.data = ndimage.sobel(data, axis=-1, mode='nearest')
    sobel_y.data = ndimage.sobel(data, axis=-2, mode='nearest')
    return create_vf([sobel_x, sobel_y], VectorField.XY)

def compress(data, upper=99.8, lower=25, subtract_low=False):
    res = derive(data, 'compress()')
    values = data.values
    th_upper = np.percentile(values, upper)
    th_lower = np.percentile(values, lower)

    res.data[values > th_upper] = th_upper
    res.data[values < th_lower] = th_lower
    if subtract_low:
        res.data -= th_lower
    return res

def mask_edges(data_array, d=3):
    res = derive(
            data_array,
            f'mask_edges(d={d})'
            )
    if isinstance(d, Number):
        dx = d
        dy = d
    else:
        dx = d[0]
        dy = d[1]
    res.data[:dy, :] = 0
    res.data[-dy:, :] = 0
    res.data[:, :dx] = 0
    res.data[:, -dx:] = 0
    return res

def subtract_noise(data_array, percentile=85):
    res = derive(
            data_array,
            f'subtract_noise(percentile={percentile})'
            )

    mask = np.percentile(res.data, percentile)
    noise = np.average(res.data[res.data < mask])
    res.data = res.data - noise
    return res

def attenuate_edges(data_array):
    ny, nx = data_array.shape
    x = np.abs(np.linspace(-1, 1, nx))
    y = np.abs(np.linspace(-1, 1, ny))
    xv, yv = np.meshgrid(x,y)
    mask = np.max([xv,yv], axis=0)
    mask = 1 - mask**4
    res = derive(
            data_array,
            f'attenuate_edges()'
            )
    res.data *= mask
    return res

