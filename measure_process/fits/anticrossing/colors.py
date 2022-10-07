import numpy as np
import xarray as xr

from projects.xarray_utils.vector_fields import create_vf, vf2polar, vf_hsv_plot, VectorField

x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)

xg, yg = np.meshgrid(x, y)

xa = xr.DataArray(xg, coords=(y,x), dims=('y','x'))
#xa.coords = {'x':x, 'y':y}

ya = xr.DataArray(yg, coords=(y,x), dims=('y','x'))
#xy.coords = {'x':x, 'y':y}

xy = create_vf([xa,ya], VectorField.XY)

r, phi = vf2polar(xy)

xy.data[0][r > 1.0] = np.nan

vf_hsv_plot(xy, colors='dark')