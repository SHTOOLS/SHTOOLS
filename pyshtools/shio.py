"""
pyshtools Spherical Harmonic I/O, Storage, and Conversion Routines.

This submodule of pyshtools defines the following functions:

Spherical harmonic I/O
----------------------
SHRead           Read spherical harmonic coefficients from an ascii-formatted
                 file.
SHReadH          Read spherical harmonic coefficients from an ascii-formatted
                 file with a header line.
SHReadError      Read spherical harmonic coefficients and associated errors
                 from an ascii-formatted file.
SHReadErrorH     Read spherical harmonic coefficients and associated errors
                 from an ascii-formatted file with a header line.
SHRead2          Read spherical harmonic coefficients from a CHAMP or GRACE-
                 like ascii-formatted file.
SHRead2Error     Read spherical harmonic coefficients and associated errors
                 from a CHAMP or GRACE-like ascii-formatted file.
SHReadJPL        Read spherical harmonic coefficients from a JPL ascii-
                 formatted file.
SHReadJPLError   Read spherical harmonic coefficients and associated errors
                 from a JPL ascii-formatted file.
SHReadGFC        Read spherical harmonic coefficients or associated errors
                 from an ICGEM GFC ascii-formatted file.

Spherical harmonic storage
--------------------------
SHCilmToCindex   Convert a three-dimensional array of complex spherical
                 harmonic coefficients to a two-dimensional indexed array.
SHCindexToCilm   Convert a two-dimensional indexed array of complex spherical
                 harmonic coefficients to a three-dimensional array.
SHCilmToVector   Convert a 3-dimensional array of real spherical harmonic
                 coefficients to a 1-dimensional ordered array.
SHVectorToCilm   Convert a 1-dimensional indexed vector of real spherical
                 harmonic coefficients to a 3-dimensional array.
YilmIndexVector  Determine the index of a 1D ordered vector of spherical
                 harmonic coefficients corresponding to I, L, and M.

Spherical harmonic conversions
------------------------------
SHrtoc           Convert real spherical harmonics to complex form.
SHctor           Convert complex spherical harmonics to real form.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import datetime as _datetime
import numpy as _np

from ._SHTOOLS import SHRead
from ._SHTOOLS import SHReadH
from ._SHTOOLS import SHReadError
from ._SHTOOLS import SHReadErrorH
from ._SHTOOLS import SHRead2
from ._SHTOOLS import SHRead2Error
from ._SHTOOLS import SHReadJPL
from ._SHTOOLS import SHReadJPLError
from ._SHTOOLS import SHCilmToCindex
from ._SHTOOLS import SHCindexToCilm
from ._SHTOOLS import SHCilmToVector
from ._SHTOOLS import SHVectorToCilm
from ._SHTOOLS import SHrtoc
from ._SHTOOLS import SHctor


def _yyyymmdd_to_year_fraction(date):
    """Convert YYYMMDD.DD date string or float to YYYY.YYY"""
    date = str(date)
    if '.' in date:
        date, residual = str(date).split('.')
        residual = float('0.' + residual)
    else:
        residual = 0.0

    date = _datetime.datetime.strptime(date, '%Y%m%d')
    date += _datetime.timedelta(days=residual)

    year = date.year
    year_start = _datetime.datetime(year=year, month=1, day=1)
    next_year_start = _datetime.datetime(year=year + 1, month=1, day=1)
    year_duration = next_year_start - year_start

    year_elapsed = date - year_start
    fraction = year_elapsed / year_duration

    return year + fraction


def _derivative(epoch, ref_epoch, trnd, periodic):
    """Return sum of the time-variable part of the coefficients

    The formula is:
    G(t)=G(t0) + trnd*(t-t0 ) +
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


def read_icgem_gfc(filename, lmax=None, usecols=(3, 4), epoch=None):
    """Read spherical harmonic coefficients from an ICGEM GFC ascii-formatted file.

    This function only reads file with the gravity field spherical
    harmonic coefficients.

    Returns
    -------
    cilm : array
        Array with the coefficients with the shape
        (2, lmax + 1, lmax + 1) for the given epoch.
    gm : float
        Standart gravitational constant of the model, in m**3/s**2
    r0 : float
        Reference radius of the model, in meters.

    Parameters
    ----------
    filename : str
        The ascii-formatted filename containing the spherical harmonic coefficients.
    usecols : sequence, optional
        Which two columns to read. Can be used to read errors of the
        coefficients instead of the coefficients themself. Can be (3, 4) or
        (5, 6) or (7, 8). Default is (3, 4), i.e. read coefficients.
    lmax : int, optional
        Maximum degree to read from file, if None then maximum degree of the
        model will be used as well as if lmax < 0 or lmax > lmax_model.
    epoch : str or float, optional
        The epoch time to calculate time-variable coefficients in YYYYMMDD.DD
        format. If None then reference epoch t0 of the model will be used.
        If format of the file is 'icgem2.0' then epoch must be specified.
    """

    if usecols not in ((3, 4), (5, 6), (7, 8)):
        raise ValueError('Wrong column number for usecols.')

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
            raise ValueError('No standart gravitational constant in the\
                    header')

        radius = float(header['radius'])

        lmax_model = int(header['max_degree'])
        if lmax is None or lmax < 0 or lmax > lmax_model:
            lmax = lmax_model

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

# ---------------------------------------------------------------------
# --- Define a Python function that replaces the Fortran
# --- equivalent that uses different indexing conventions.
# ---------------------------------------------------------------------


def YilmIndexVector(i, l, m):
    return l**2 + (i - 1) * l + m
