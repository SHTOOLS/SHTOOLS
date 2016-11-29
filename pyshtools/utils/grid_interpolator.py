#!/usr/bin/env python

"""
2D Interpolator for mapping purposes.
"""

from scipy.ndimage.interpolation import map_coordinates as _map_coordinates
import numpy as _np


def interpolate_2dgrid(x_new, y_new, x1d_old, y1d_old, grid_old, order=2):
    # save shape of input array
    inshape = x_new.shape
    nx1d = len(x1d_old)
    ny1d = len(y1d_old)

    # add north and south pole
    shape_with_poles = grid_old.shape[0] + 2, grid_old.shape[1]
    grid_old_with_poles = _np.empty(shape_with_poles)
    grid_old_with_poles[1:-1, :] = grid_old
    grid_old[0, :] = _np.sum(grid_old[0])
    grid_old[-1, :] = _np.sum(grid_old[-1])

    # add wrapped longitude to x coordinates
    x1d_old_wrap = _np.concatenate([[x1d_old[-1] - 360.], x1d_old,
                                    [x1d_old[0] + 360.]])
    ix1d_old_wrap = _np.concatenate([[nx1d - 1], _np.arange(nx1d),
                                     [0]])

    # find index coordinate of flattened x_new, y_new
    ix = _np.interp(x_new.flat, x1d_old_wrap, ix1d_old_wrap)
    iy = _np.interp(y_new.flat, y1d_old, _np.arange(ny1d))

    # make index grid
    points = _np.array([ix, iy])

    # get flattened values from grid
    values = _map_coordinates(grid_old, points, order=order, mode='wrap')

    # reshape to original shape
    values = values.reshape(inshape)
    return values
