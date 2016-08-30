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
    other - Other routines.

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

# ---- Import all wrapped SHTOOLS functions into shtools submodule
from . import _SHTOOLS as shtools

# ---- Import classes into pyshtools namespace
from . import shclasses
from .shclasses import SHCoeffs, SHGrid, SHWindow

# ---- Import shtools submodules ----
from . import constant
from . import legendre
from . import expand
from . import shio
from . import spectralanalysis
from . import localizedspectralanalysis
from . import rotate
from . import gravmag
from . import other


# ---- Bind two new functions to the list of all shtools routines ----
_SHTOOLS.PlmIndex = legendre.PlmIndex
_SHTOOLS.YilmIndexVector = shio.YilmIndexVector


# ---------------------------------------------------------------------
# ---- Fill the pyshtools module doc strings and pyshtools constant
# ---- infostrings with documentation from external files. The doc files
# ---- are generated during intitial compilation of pyshtools from md
# ---- formatted text files.
# ---------------------------------------------------------------------
_pydocfolder = _os.path.abspath(_os.path.join(_os.path.dirname(__file__),
                                              'doc'))

for _name, _func in _SHTOOLS.__dict__.items():
    if callable(_func):
        try:
            _path = _os.path.join(_pydocfolder, _name.lower() + '.doc')

            with open(_path) as _pydocfile:
                _pydoc = _pydocfile.read()

            _func.__doc__ = _pydoc
        except IOError as msg:
            print(msg)

for _name in _constant.planetsconstants.__dict__.keys():
    try:
        _path = _os.path.join(_pydocfolder, 'constant_' + _name.lower() +
                              '.doc')

        with open(_path) as _pydocfile:
            _pydoc = _pydocfile.read()

        setattr(getattr(constant, _name), '_infostring', _pydoc)

    except IOError as msg:
        print(msg)


# ---- Rewrite shtools doc string ----
shtools.__doc__ = (
    'pyshtools submodule that includes all pyshtools routines, with the\n' +
    'exception of shclasses.')


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['constant', 'shclasses', 'SHCoeffs', 'SHGrid', 'SHWindow',
           'legendre', 'expand', 'shio', 'spectralanalysis',
           'localizedspectralanalysis', 'rotate', 'gravmag', 'other']
