"""
pyshtools
=========

pyshtools a scientific package that can be used to perform spherical harmonic
transforms and reconstructions, rotations of data expressed in spherical
harmonics, and multitaper spectral analyses on the sphere.

This module imports the following classes and subpackages into the
main namespace:

    SHCoeffs          : Class for spherical harmonic coefficients.
    SHGrid            : Class for global grids.
    SHWindow          : Class for localized spectral analyses.
    Slepian           : Class for Slepian functions.
    SHGravCoeffs      : Class for gravitational potential spherical harmonic
                        coefficients.
    SHMagCoeffs       : Class for magnetic potential spherical harmonic
                        coefficients.

    shclasses         : All pyshtools classes and subclasses.
    shtools           : All Python-wrapped Fortran 95 routines.
    constants         : pyshtools constants.
    legendre          : Legendre functions.
    expand            : Spherical harmonic expansion routines.
    shio              : Spherical harmonic I/O, storage, and conversion
                        routines.
    spectralanalysis  : Global and localized spectral analysis routines.
    rotate            : Spherical harmonic rotation routines.
    gravmag           : Gravity and magnetics routines.
    utils             : pyshtools utilities.
    backends          : Routines for selecting the spherical harmonic
                        transform backend.

For further information, consult the web documentation at

   https://shtools.github.io/SHTOOLS/

and the GitHub project page at

   https://github.com/SHTOOLS/SHTOOLS
"""
from ._version import get_versions as _get_versions

# ---- Import shtools subpackages ----
from . import backends
from . import constants
from . import datasets
from . import shclasses
from . import legendre
from . import expand
from . import shio
from . import spectralanalysis
from . import rotate
from . import gravmag
from . import utils
from .backends import shtools

# ---- Import principal classes into pyshtools namespace
from .shclasses import SHCoeffs
from .shclasses import SHGrid
from .shclasses import SHWindow
from .shclasses import Slepian
from .shclasses import SHGravCoeffs
from .shclasses import SHMagCoeffs

__version__ = _get_versions()["version"]
__commit__ = _get_versions()["full-revisionid"]
__author__ = 'SHTOOLS developers'

# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['constants', 'shclasses', 'legendre', 'expand', 'shio', 'shtools',
           'spectralanalysis', 'rotate', 'gravmag', 'utils', 'backends',
           'SHCoeffs', 'SHGrid', 'SHWindow', 'Slepian', 'SHGravCoeffs',
           'SHMagCoeffs', 'datasets']
