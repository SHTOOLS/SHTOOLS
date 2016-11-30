"""
ICGEM-format read support
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division

import numpy as _np

from pyshtools.utils.datetime import _yyyymmdd_to_year_fraction


def _time_variable_part(epoch, ref_epoch, trnd, periodic):
    """Return sum of the time-variable part of the coefficients

    The formula is:
    G(t) = G(t0) + trnd*(t-t0) +
        asin1*sin(2pi/p1 * (t-t0)) + acos1*cos(2pi/p1 * (t-t0)) +
        asin2*sin(2pi/p2 * (t-t0)) + acos2*cos(2pi/p2 * (t-t0))

    This function computes all terms after G(t0).
    """
    delta_t = epoch - ref_epoch
    trend = trnd * delta_t
    periodic_sum = _np.zeros_like(trnd)
    for period in periodic:
        for trifunc in periodic[period]:
            coeffs = periodic[period][trifunc]
            if trifunc == 'acos':
                periodic_sum += coeffs * _np.cos(2 * _np.pi / period * delta_t)
            elif trifunc == 'asin':
                periodic_sum += coeffs * _np.sin(2 * _np.pi / period * delta_t)
    return trend + periodic_sum


def read_icgem_gfc(filename, errors=None, lmax=None, epoch=None):
    """Read spherical harmonic coefficients from an ICGEM GFC ascii-formatted file.

    This function only reads files with the gravity field spherical
    harmonic coefficients.

    Returns
    -------
    cilm : array
        Array with the coefficients with the shape (2, lmax + 1, lmax + 1)
        for the given epoch.
    gm : float
        Standard gravitational constant of the model, in m**3/s**2.
    r0 : float
        Reference radius of the model, in meters.
    errors : array, optional
        Array with the errors of the coefficients with the shape
        (2, lmax + 1, lmax + 1) for the given epoch.

    Parameters
    ----------
    filename : str
        The ascii-formatted filename containing the spherical harmonic coefficients.
    errors : str, optional
        Which errors to read. Can be either "calibrated", "formal" or
        None. Default is None.
    lmax : int, optional
        Maximum degree to read from the file. If lmax is None, less than 0, or
        greater than lmax_model, the maximum degree of the model will be used.
    epoch : str or float, optional
        The epoch time to calculate time-variable coefficients in YYYYMMDD.DD
        format. If None then reference epoch t0 of the model will be used.
        If format of the file is 'icgem2.0' then epoch must be specified.
    """

    # read header
    header = {}
    header_keys = ['modelname', 'product_type', 'earth_gravity_constant',
                   'gravity_constant', 'radius', 'max_degree', 'errors',
                   'tide_system', 'norm', 'format']

    with open(filename, 'r') as f:
        for line in f:
            if 'end_of_head' in line:
                break
            for key in header_keys:
                if key in line:
                    header[key] = line.strip().split()[1]

        if header['product_type'] != 'gravity_field':
            raise ValueError('This function reads only gravity_field data product.')

        is_v2 = False
        if 'format' in header and header['format'] == 'icgem2.0':
            is_v2 = True

        if epoch is None and is_v2:
            raise ValueError('Epoch must be specified for the "icgem2.0" format.')
        elif epoch is not None:
            epoch = _yyyymmdd_to_year_fraction(epoch)

        if 'earth_gravity_constant' in header:
            gravity_constant = float(header['earth_gravity_constant'])
        elif 'gravity_constant' in header:
            gravity_constant = float(header['gravity_constant'])
        else:
            raise ValueError('No standard gravitational constant in the header.')

        radius = float(header['radius'])

        lmax_model = int(header['max_degree'])
        if lmax is None or lmax < 0 or lmax > lmax_model:
            lmax = lmax_model

        if errors is not None:
            valid_err = ('calibrated', 'formal', 'calibrated_and_formal')
            if header['errors'] == 'no':
                raise ValueError('This model has no errors.')
            elif errors not in valid_err[:-1]:
                raise ValueError('Errors can be either "formal", "calibrated" or None.')
            elif header['errors'] in valid_err and errors in valid_err[:-1]:
                if (errors, header['errors']) == valid_err[1:]:
                    err_cols = (7, 8)
                elif header['errors'] != errors:
                    raise ValueError('This model has no {} errors.'.format(errors))
                else:
                    err_cols = (5, 6)

        cilm = _np.tile(_np.zeros((lmax + 1, lmax + 1)), (4, 1, 1))
        ref_epoch = _np.zeros((lmax + 1, lmax + 1))
        trnd = _np.zeros_like(cilm)
        periodic = {}

        # read coefficients
        for line in f:
            line = line.replace('D', 'E').strip().split()

            l, m = int(line[1]), int(line[2])
            if m > lmax:
                break
            if l > lmax:
                continue

            key = line[0]

            value_cs = [float(line[3]), float(line[4]), 0, 0]
            if errors:
                value_cs[2:] = float(line[err_cols[0]]),\
                    float(line[err_cols[1]])

            if key == 'gfc':
                cilm[:, l, m] = value_cs
            elif key == 'gfct':
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-2])
                    t1i = _yyyymmdd_to_year_fraction(line[-1])
                    if not t0i <= epoch < t1i:
                        continue
                else:
                    t0i = _yyyymmdd_to_year_fraction(line[-1])

                cilm[:, l, m] = value_cs
                ref_epoch[l, m] = t0i
            elif key == 'trnd':
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-2])
                    t1i = _yyyymmdd_to_year_fraction(line[-1])
                    if not t0i <= epoch < t1i:
                        continue
                trnd[:, l, m] = value_cs
            elif key in ('acos', 'asin'):
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-3])
                    t1i = _yyyymmdd_to_year_fraction(line[-2])
                    if not t0i <= epoch < t1i:
                        continue

                period = float(line[-1])
                if period not in periodic:
                    arr = _np.zeros_like(cilm)
                    periodic[period] = {'acos': arr,
                                        'asin': arr.copy()}

                periodic[period][key][:, l, m] = value_cs

    if epoch is None:
        epoch = ref_epoch

    cilm += _time_variable_part(epoch, ref_epoch, trnd, periodic)

    if errors:
        return cilm[:2], gravity_constant, radius, cilm[2:]
    else:
        return cilm[:2], gravity_constant, radius
