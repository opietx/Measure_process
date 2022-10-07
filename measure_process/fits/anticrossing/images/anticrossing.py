import numpy as np
import scipy.optimize as optimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as pt

from projects.xarray_utils.data import derive, xy_coord_values

from .fit_utils import ModelParameters, ModelParameter
from .filters import mask_edges
import numpy as np
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def xy2ixy(data, xv, yv):
    x,y = xy_coord_values(data)
    nx,ny = len(x), len(y)
    # ix = interp1d(x, np.arange(nx))(xv)
    ix = find_nearest(xv,data.P6.data)
    iy = find_nearest(xv,data.P5.data)
    # iy = interp1d(y, np.arange(ny))(yv)
    return ix, iy


def ixy2xy(data, ix, iy):
    x,y = xy_coord_values(data)
    nx,ny = len(x), len(y)
    xv = interp1d(np.arange(nx), x, fill_value="extrapolate")(ix)
    yv = interp1d(np.arange(ny), y, fill_value="extrapolate")(iy)
    return xv, yv


def pt_line(p1, p2, color):
    pt.plot([p1[0],p2[0]], [p1[1],p2[1]], '.-', color=color)


def plot_cross(points):
    p1,p2,p3,p4,p5,p6 = list(points)

    pt.plot((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, 'ko')
    pt_line(p1, p2, 'lime')
    pt_line(p1, p3, 'orange')
    pt_line(p1, p4, 'orange')
    pt_line(p2, p5, 'orange')
    pt_line(p2, p6, 'orange')


def add_line(img, p1, p2, sigma, fade=True, sign=1, vf=False):
    p1 = np.array(p1)
    p2 = np.array(p2)
    pp = p2 - p1
    norm_pp = np.linalg.norm(pp)
    pp = pp/norm_pp
    g1 = img - p1
    # distance to infinite line through p1-p2
    gd = np.cross(g1, pp)
    # relative point on line p1 - p2: 0.0 is next to p1, 1.0 is next to p2
    # explicitly writting the dot product is much faster than np.dot..
#    s_pp = np.dot(g1, pp) / norm_pp # this is slow...
    s_pp = (g1[:,:,0] * pp[0] + g1[:,:,1] * pp[1]) / norm_pp

    # Gaussian shape
    if fade:
        # fade towards the end with cos shape: smooth 1 at stort, smooth 0 at end
        g = ((s_pp > 0) & (s_pp < 1)) * np.exp(gd**2*(-0.5/sigma)) * (0.5+0.5*np.cos(s_pp*np.pi))
    else:
        g = ((s_pp > 0) & (s_pp < 1)) * np.exp(gd**2*(-0.5/sigma))

    if vf:
        # dx,dy is direction perpendicular to pp.
        # This is direction of the gradient in CSD
        dx,dy = -pp[1],pp[0]
        return np.array([dx*g, dy*g])*sign
    else:
        return g


def make_cross(x, y, sigma, legs, leg_length, vf,
               center_x, center_y, width, angle_left, angle_right, angle_center):
    nx = len(x)
    ny = len(y)

    center = np.array([center_x, center_y])
    d = width * np.array([np.cos(angle_center), np.sin(angle_center)])
    lx = nx*leg_length
    ly = ny*leg_length
    ar = angle_center + angle_right
    al = angle_center + angle_left
    dr = np.array([lx*np.cos(ar), ly*np.sin(ar)])
    dl = np.array([lx*np.cos(al), ly*np.sin(al)])

    p1 = center + d/2
    p2 = center - d/2
    p3 = p1 + dr
    p4 = p1 + dl
    p5 = p2 - dr
    p6 = p2 - dl

    g = np.zeros((ny,nx,2))
    g[:,:,0] = (np.arange(nx)-(nx/2))[None,:]
    g[:,:,1] = (np.arange(ny)-(ny/2))[:,None]

    lines = []
    if legs is None or 'center' in legs or 'c' in legs:
        lines.append( add_line(g, p1, p2, sigma, fade=False, vf=vf) )
    if legs is None or 'up' in legs or 2 in legs:
        lines.append( add_line(g, p1, p3, sigma, sign=1, vf=vf) )
    if legs is None or 'up' in legs or 4 in legs:
        lines.append( add_line(g, p2, p5, sigma, sign=-1, vf=vf) )
    if legs is None or 'down' in legs or 1 in legs:
        lines.append( add_line(g, p1, p4, sigma, sign=-1, vf=vf) )
    if legs is None or 'down' in legs or 3 in legs:
        lines.append( add_line(g, p2, p6, sigma, sign=1, vf=vf) )

    if vf:
        img = np.sum(lines, axis=(0))
    else:
        # use maximum of gradient magnitudes
        img = np.max(lines, axis=0)

#    with np.printoptions(precision=3):
#        print(np.array([p1,p2,p3,p4,p5,p6]))

    return img, [p1,p2,p3,p4,p5,p6]


def match_img(img, x, y, sigma_cross, legs, leg_length, vf, model_values):
    img_cross, points = make_cross(x, y, sigma_cross, legs, leg_length, vf, **model_values)
    norm_img = np.linalg.norm(img)
    norm_cross = np.linalg.norm(img_cross)
    match = np.dot(img_cross.flat, img.flat)/(norm_img*norm_cross)
    return match


def cross_cost_func(fit_values, model_params, x, y, img, sigma, legs, leg_length, vf, fit_xy_rotation):
    param_values = model_params.get_model_values(fit_values)

    model_values = param_values.copy()
    center_x, center_y = rot2D(-fit_xy_rotation) @ [model_values['center_x'], model_values['center_y']]
    model_values['center_x'] = center_x
    model_values['center_y'] = center_y

    match = match_img(img, x, y, sigma, legs, leg_length, vf, model_values)

    # add penalty for difference between angles.
    # angle should be approximately symmetric around 0
    angle_right = model_values['angle_right']
    angle_left = model_values['angle_left']
    angle_diff = abs(angle_right + angle_left)/np.pi
    if vf:
        angle_penalty = (1- 0.3 * angle_diff)
    else:
        angle_penalty = (1- 0.4 * angle_diff)
    return -match * angle_penalty


def rot2D(angle):
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array([[c,s],[-s,c]])


def fit_anti_crossing(xr_data, fit_param_names, sigma, legs, leg_length, start_x, start_y, vf,
                      fit_xy_rotation=0, center_angle=0.5*np.pi):
    '''
    Fits the anti-crossing using non-linear optimization fit using a model
    of the anti-crossing with addition lines.
    Optimization is done with Sequential Least Squares Procedure.

    Args:
        xr_data: DataArray with gradient of CSD scan
        fit_param_names: list with strings of parameters to fit.
        sigma: standard deviation of Gaussian curve used in the model
        legs: legs to use for model: ['c', 1, 2, 3, 4]
        leg_length: length of the legs of the model in fraction of the image size
        start_x: initial value for x-coordinate of anti-crossing center
        start_y: initial value for y-coordinate of anti-crossing center
        vf: if True uses model with x and y component of gradient (vector field)
        fit_xy_rotation: angle for rotation x,y coordinate in model. This might improve fit quality.
        center_angle: angle of the anti-crossing line

    Returns:
        DataArray with fitted model and results in attrs dictionary.
    '''
    x,y = xy_coord_values(xr_data)
    img = xr_data.data
    nx = len(x)
    ny = len(y)


    # image coordinates (-ny/2:+ny/2, -nx/2:+nx/2)

    # convert start_x, start_y to image coordinates
    # print(xr_data)
    start_ix, start_iy = xy2ixy(xr_data, start_x, start_y)
    # print('hereeee', start_ix, start_iy)
    start_ix -= nx/2
    start_iy -= ny/2
    if not (-0.4*nx < start_ix < 0.4*nx and -0.4*ny < start_iy < 0.4*ny):
        print(f'Bad start condition {start_x:5.2f}, {start_y:5.2f}')
        start_ix, start_iy = 0.0, 0.0

    init_x, init_y = rot2D(fit_xy_rotation) @ [start_ix, start_iy]

    # parameters in image coordinates
    param_dict = {
        'center_x': ModelParameter(init_x, -0.4*nx, 0.4*nx),
        'center_y': ModelParameter(init_y, -0.4*ny, 0.4*ny),
        'width': ModelParameter(max(3.0, 0.05*ny), 0.5, 0.7*ny),
        'angle_center': ModelParameter(center_angle,
                                       center_angle - 0.05*np.pi,
                                       center_angle + 0.05*np.pi),
        'angle_right': ModelParameter(-np.pi*0.25, -np.pi*0.35, -np.pi*0.15),
        'angle_left': ModelParameter(np.pi*0.25, np.pi*0.15, np.pi*0.35),
        }

    model_params = ModelParameters(param_dict, fit_param_names)

    res = optimize.minimize(
            cross_cost_func,
            model_params.get_initial_values(),
            args=(model_params, x, y, img, sigma, legs, leg_length, vf, fit_xy_rotation),
            method='SLSQP',
            options={
                'maxiter':100,
                'eps': 1e-4,
                'ftol': 1e-7,
                },
            bounds=model_params.get_boundaries()
            )

#    print(res.nit, res.nfev)
#    print(np.array2string(res.jac, precision=3))

    param_values = model_params.get_model_values(res.x)
    model_values = param_values.copy()
    center_x, center_y = rot2D(-fit_xy_rotation) @ [model_values['center_x'], model_values['center_y']]
    model_values['center_x'] = center_x
    model_values['center_y'] = center_y

    img_cross, points = make_cross(x, y, sigma, legs, leg_length, vf, **model_values)
    match = match_img(img, x, y, sigma, ['c',1,2,3,4], leg_length, vf, model_values)
    match_legs = {}
    for leg in [1,2,3,4]:
        match_legs[leg] = match_img(img, x, y, sigma, [leg], leg_length, vf, model_values)

    # convert image coordinates to scan coordinates
    points_xy = []
    for point in points:
        points_xy.append(np.array(ixy2xy(xr_data, point[0]+nx/2, point[1]+ny/2)))

    c_ix = model_values['center_x']
    c_iy = model_values['center_y']
    cx, cy = ixy2xy(xr_data, c_ix+nx/2, c_iy+ny/2)

    result = derive(xr_data, 'AntiCrossingFit()')
    result.data = img_cross
    attrs = result.attrs
    attrs['success'] = res.success
    attrs['match'] = match
    attrs['match_legs'] = match_legs
    attrs['points'] = np.array(points_xy)
    attrs['center'] = np.array([cx, cy])
    attrs['width'] = np.linalg.norm(points_xy[1] - points_xy[0])
    attrs['fitted_model'] = model_values

#    model_values = model_params.get_model_values(model_params.get_initial_values())
#    cross_init, points_init = make_cross(x, y, sigma, legs, leg_length, vf, **model_values)
#    points_xy = []
#    for point in points_init:
#        points_xy.append(np.array(ixy2xy(xr_data, point[0]+nx/2, point[1]+ny/2)))
#    points_init = points_xy
#    attrs['cross_init'] = cross_init
#    attrs['points_init'] = points_init

    return result


def estimate_center(img_mag, img_angle, grad_angle=0.0, both_dir=False):
    '''
    Estimates the center of the anti-crossing.
    Searching the position of the maximum magnitude of the gradient with
    the specified direction. A margin of 0.15*PI is applied to the direction.

    Args:
        img_amp: 2D xarray DataArray with magnitude of gradient.
        img_angle: 2D xarray DataArray with direction of gradient.
        grad_angle: direction of the gradient.
        both_dir: if True the opposite angle is also accepted.

    Returns:
        x-coordinate of center
        y-coordinate of center
        image with magnitude of
    '''
    x,y = xy_coord_values(img_mag)
    nx, ny = len(x), len(y)
    th = np.percentile(img_mag, 70)
    img_ac = derive(img_mag)

    angle = img_angle.values
    # rotate setting searched angle to 0.0 in units of PI
    angle = (img_angle.values - grad_angle) / np.pi
    # normalize angles to +/- 1, or +/- 1/2 for both directions
    if both_dir:
        angle = (angle + 2.5) % 1 - 0.5
    else:
        angle = (angle + 3) % 2 - 1

    # keep only gradient with direction of anti-crossing and above threshold
    img_ac.data[(img_mag.values < th)
                | (angle < -0.2)
                | (angle > +0.2)] = 0.0

    # search maximimum value excluding margin of 20% on each side
    img_ac = mask_edges(img_ac, d=(20*nx//100, 20*ny//100))
    index_ac = np.unravel_index(np.argmax(img_ac.data), img_mag.shape)
    ac_x = x[index_ac[1]]
    ac_y = y[index_ac[0]]
    return ac_x, ac_y, img_ac
