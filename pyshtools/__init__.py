"""
pyshtools
=========

pyshtools is an archive of scientific routines that can be used to
perform spherical harmonic transforms and reconstructions, rotations
of data expressed in spherical harmonics, and multitaper spectral
analyses on the sphere.

This module imports the following classes and subpackages into the
main namespace:

    SHCoeffs - A high level class for spherical harmonic coefficients.
    SHGrid - A high level classes for global grids.
    SHWindow - A high level classes for localization windows.
    shclasses - All pyshtools classes and subclasses.
    shtools - All Python-wrapped Fortran 95 routines.
    constant - pyshtools constants.
    legendre - Legendre functions.
    expand - Spherical harmonic expansion routines.
    shio - Spherical harmonic I/O, storage, and conversion routines.
    spectralanalysis - Global and localized spectral analysis routines.
    rotate - Spherical harmonic rotation routines.
    gravmag - Gravity and magnetics routines.
    utils - pyshtools utilities.

For further information, consult the web documentation at

   https://shtools.oca.eu/

and the GitHub project page at

   https://github.com/SHTOOLS/SHTOOLS
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

__version__ = '4.2'
__author__ = 'SHTOOLS developers'

import os as _os
import numpy as _np

# ---- Import shtools subpackages ----
from . import shtools
from . import constant
from . import shclasses
from . import legendre
from . import expand
from . import shio
from . import spectralanalysis
from . import rotate
from . import gravmag
from . import utils

# ---- Import classes into pyshtools namespace
from .shclasses import SHCoeffs, SHGrid, SHWindow

# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['constant', 'shclasses', 'SHCoeffs', 'SHGrid', 'SHWindow',
           'legendre', 'expand', 'shio', 'spectralanalysis',
           'rotate', 'gravmag', 'utils']
