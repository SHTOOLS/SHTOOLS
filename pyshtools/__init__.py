"""
pyshtools
=========

pyshtools is an archive of scientific routines that can be used to
perform spherical harmonic transforms and reconstructions, rotations
of data expressed in spherical harmonics, and multitaper spectral
analyses on the sphere.

This module imports the following classes and submodules into the
main namespace:

    SHCoeffs - A high level class for spherical harmonic coefficients.
    SHGrid - A high level classes for global grids.
    SHWindow - A high level classes for localization windows.
    shclasses - All pyshtools classes and subclasses.
    shtools - All pyshtools routines.
    legendre - Legendre functions.
    expand - Spherical harmonic expansion routines.
    shio - Spherical harmonic I/O, storage, and conversion routines.
    spectralanalysis - Global spectral analysis routines.
    localizedspectralanalysis - Localized spectral analysis routines.
    rotate - Spherical harmonic rotation routines.
    gravmag - Gravity and magnetics routines.
    constant - pyshtools constants.
    utils - Other routines.

For further information, consult the web documentation at

   http://shtools.ipgp.fr/

and the GitHub project page at

   https://github.com/SHTOOLS/SHTOOLS
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

__version__ = '3.4'
__author__ = 'SHTOOLS developers'

import os as _os
import numpy as _np


# ---- Import classes into pyshtools namespace
from . import shclasses
from .shclasses import SHCoeffs, SHGrid, SHWindow

# ---- Import shtools subpackages ----
from . import shtools
from . import constant
from . import legendre
from . import expand
from . import shio
from . import utils

# ---- Import shtools submodules ----
from . import spectralanalysis
from . import localizedspectralanalysis
from . import rotate
from . import gravmag


# ---- Check the exit status of Fortran routines, raise exceptions, and 
# ---- strip exitstatus from the Python return values.
class SHToolsError(Exception):
    pass


def _shtools_status_message(status):
    '''
    Determine error message to print when a SHTOOLS Fortran 95 routine exits
    improperly.
    '''
    if (status == 1):
        errmsg = 'Improper dimensions of input array.'
    elif (status == 2):
        errmsg = 'Improper bounds for input variable.'
    elif (status == 3):
        errmsg = 'Error allocating memory.'
    elif (status == 4):
        errmsg = 'File IO error.'
    else:
        errmsg = 'Unhandled Fortran 95 error.'
    return errmsg


def _raise_errors(func):
    def wrapped_func(*args, **kwargs):
        returned_values = func(*args, **kwargs)
        if returned_values[0] != 0:
            raise SHToolsError(_shtools_status_message(returned_values[0]))
        elif len(returned_values) == 2:
            return returned_values[1]
        else:
            return returned_values[1:]
    wrapped_func.__doc__ = func.__doc__
    return wrapped_func


_fortran_subroutines = (legendre._fortran_subroutines +
                        expand._fortran_subroutines)


for _func in _fortran_subroutines:
    setattr(_SHTOOLS, _func, _raise_errors(getattr(_SHTOOLS, _func)))


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['constant', 'shclasses', 'SHCoeffs', 'SHGrid', 'SHWindow',
           'legendre', 'expand', 'shio', 'spectralanalysis',
           'localizedspectralanalysis', 'rotate', 'gravmag', 'utils']
