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

    # add north/south pole and wrapped longitude values to the grid
    shape_wrapped = grid_old.shape[0] + 2, grid_old.shape[1] + 2
    grid_wrapped = _np.empty(shape_wrapped)
    grid_wrapped[1:-1, 1:-1] = grid_old
    grid_wrapped[0, :] = _np.mean(grid_old[0, :])
    grid_wrapped[-1, :] = _np.mean(grid_old[-1, :])
    grid_wrapped[1:-1, 0] = grid_old[:, -1]
    grid_wrapped[1:-1, -1] = grid_old[:, 0]

    # add extra poles and wrapped longitudes coordinate arrays
    x1d_wrapped = _np.concatenate([[-90.], x1d_old, [90.]])
    ix1d_wrapped = _np.arange(nx1d + 2)

    y1d_wrapped = _np.concatenate([[y1d_old[-1] - 360.], y1d_old,
                                    [y1d_old[0] + 360.]])
    iy1d_wrapped = _np.arange(ny1d + 2)

    # find index coordinate of flattened x_new, y_new
    ix = _np.interp(x_new.flat, x1d_wrapped, ix1d_wrapped)
    iy = _np.interp(y_new.flat, y1d_wrapped, iy1d_wrapped)

    # make index grid
    points = _np.array([ix, iy])

    # get flattened values from grid
    values = _map_coordinates(grid_wrapped, points, order=order,
                              mode='wrap')

    # reshape to original shape
    values = values.reshape(inshape)
    return values
