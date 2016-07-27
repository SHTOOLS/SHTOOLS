"""
pyshtools defines several classes that facilitate the interactive
examination of geographical gridded data and their associated
spherical harmonic coefficients. Subclasses are used to handle different
internal data types and superclasses are used to implement interface
functions and documentation.

pyshtools class structure:

    SHCoeffs
        SHRealCoeffs
        SHComplexCoeffs

    SHGrid
        DHRealGrid
        DHComplexGrid
        GLQRealGrid
        GLQComplexGrid

    SHWindow
        SHSymmetricWindow
        SHAsymmetricWindow

For more information, see the documentation for the top level classes.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

from .. import _SHTOOLS as _shtools

from .SHCoeffs import SHCoeffs, SHRealCoeffs, SHComplexCoeffs
from .SHGrid import (SHGrid, DHRealGrid, DHComplexGrid, GLQRealGrid,
                     GLQComplexGrid)
from .SHWindow import SHWindow, SHSymmetricWindow, SHAsymmetricWindow
