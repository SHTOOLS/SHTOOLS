"""
ICGEM-format read support
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division

import numpy as _np

from pyshtools.shio.common import _derivative
from pyshtools.utils.datetime import _yyyymmdd_to_year_fraction

def read_icgem_gfc(filename, lmax=None, errors=None, epoch=None):
    """Read spherical harmonic coefficients from an ICGEM GFC ascii-formatted file.

    This function only reads file with the gravity field spherical
    harmonic coefficients.

    Returns
    -------
    cilm : array
        Array with the coefficients with the shape
        (2, lmax + 1, lmax + 1) for the given epoch. Returns errors
        coefficients if errors parameter is set.
    gm : float
        Standart gravitational constant of the model, in m**3/s**2
    r0 : float
        Reference radius of the model, in meters.

    Parameters
    ----------
    filename : str
        The ascii-formatted filename containing the spherical harmonic coefficients.
    errors : str, optional
        Which errors to read. Can be either "calibrated", "formal" or
        None. Default is None.
    lmax : int, optional
        Maximum degree to read from file, if None then maximum degree of the
        model will be used as well as if lmax < 0 or lmax > lmax_model.
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
            raise ValueError('No standart gravitational constant in the header')

        radius = float(header['radius'])

        lmax_model = int(header['max_degree'])
        if lmax is None or lmax < 0 or lmax > lmax_model:
            lmax = lmax_model

        valid_err = ('calibrated', 'formal', 'calibrated_and_formal')
        if header['errors'] == 'no' and errors is not None:
            raise ValueError('This model has no errors.')
        elif errors not in valid_err[:-1] + (None,):
            raise ValueError('Errors can be either "formal", "calibrated" or None')
        elif header['errors'] in valid_err and errors in valid_err[:-1]:
            if (errors, header['errors']) == valid_err[1:]:
                usecols = (7, 8)
            elif header['errors'] != errors:
                raise ValueError('This model has no {} errors.'.format(errors))
            else:
                usecols = (5, 6)
        else:
            usecols = (3, 4)

        cilm = _np.tile(_np.zeros((lmax + 1, lmax + 1)), (2, 1, 1))
        ref_epoch = _np.zeros((lmax + 1, lmax + 1))
        trnd = _np.zeros_like(cilm)
        periodic = {}

        # read coefficients
        for line in f:
            line = line.replace('D', 'E').strip().split()
            key = line[0]

            l, m = int(line[1]), int(line[2])
            if l > lmax:
                break

            value_c, value_s = float(line[usecols[0]]),\
                float(line[usecols[1]])

            if key == 'gfc':
                cilm[0, l, m], cilm[1, l, m] = value_c, value_s
            elif key == 'gfct':
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-2])
                    t1i = _yyyymmdd_to_year_fraction(line[-1])
                    if not t0i <= epoch < t1i:
                        continue
                else:
                    t0i = _yyyymmdd_to_year_fraction(line[-1])

                cilm[0, l, m], cilm[1, l, m] = value_c, value_s
                ref_epoch[l, m] = t0i
            elif key == 'trnd':
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-2])
                    t1i = _yyyymmdd_to_year_fraction(line[-1])
                    if not t0i <= epoch < t1i:
                        continue
                trnd[0, l, m], trnd[1, l, m] = value_c, value_s
            elif key in ('acos', 'asin'):
                period = float(line[-1])
                if period not in periodic:
                    arr = _np.zeros_like(cilm)
                    periodic[period] = {'acos': arr,
                                        'asin': arr.copy()}
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-3])
                    t1i = _yyyymmdd_to_year_fraction(line[-2])
                    if not t0i <= epoch < t1i:
                        continue
                periodic[period][line[0]][0, l, m] = value_c
                periodic[period][line[0]][1, l, m] = value_s
            else:
                continue

    if epoch is None:
        epoch = ref_epoch

    cilm += _derivative(epoch, ref_epoch, trnd, periodic)

    return cilm, gravity_constant, radius

