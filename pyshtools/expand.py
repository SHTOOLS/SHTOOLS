"""
pyshtools Spherical Harmonic Expansion Routines.

This submodule of pyshtools defines the following functions:

Equally sampled (N by N) and equally spaced (N by 2N) Grids
-----------------------------------------------------
SHExpandDH    Expand an equally sampled or equally spaced map into spherical
              harmonics using Driscoll and Healy's (1994) sampling theorem.
MakeGridDH    Create a 2D map from a set of spherical harmonic coefficients
              that conforms with Driscoll and Healy's (1994) sampling theorem.
SHExpandDHC   Expand an equally sampled or equally spaced complex map into
              complex spherical harmonics using Driscoll and Healy's (1994)
              sampling theorem.
MakeGridDHC   Create a 2D complex map from a set of complex spherical harmonic
              coefficients that conforms with Driscoll and Healy's (1994)
              sampling theorem.

Gauss-Legendre quadrature grids
-------------------------------
SHGLQ         Precompute the weights and nodes used in the GLQ-based spherical
              harmonics routines.
SHExpandGLQ   Expand a 2D map sampled on the Gauss-Legendre quadrature nodes
              into spherical harmonics.
MakeGridGLQ   Create a 2D map from a set of spherical harmonic coefficients
              sampled on a the Gauss-Legendre quadrature nodes.
SHExpandGLQC  Expand a 2D complex map sampled on the Gauss-Legendre quadrature
              nodes into complex spherical harmonics.
MakeGridGLQC  Create a 2D complex map from a set of complex spherical harmonic
              coefficients sampled on a the Gauss-Legendre quadrature nodes.
GLQGridCoord  Compute the latitude and longitude coordinates used in Gauss-
              Legendre quadrature grids.

Other
-----
SHExpandLSQ   Expand a set of irregularly sampled data points into spherical
              harmonics using a least squares inversion.
MakeGrid2D    Create a 2D cylindrical map with arbitrary grid spacing from a
              set of spherical harmonic coefficients.
MakeGridPoint Evaluate a function expressed in spherical harmonics at a
              single point.
SHMultiply    Multiply two functions and determine the spherical harmonic
              coefficients of the resulting function.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

from ._SHTOOLS import SHExpandDH
from ._SHTOOLS import MakeGridDH
from ._SHTOOLS import SHExpandDHC
from ._SHTOOLS import MakeGridDHC
from ._SHTOOLS import SHGLQ
from ._SHTOOLS import SHExpandGLQ
from ._SHTOOLS import MakeGridGLQ
from ._SHTOOLS import SHExpandGLQC
from ._SHTOOLS import MakeGridGLQC
from ._SHTOOLS import GLQGridCoord
from ._SHTOOLS import SHExpandLSQ
from ._SHTOOLS import MakeGrid2D
from ._SHTOOLS import MakeGridPoint
from ._SHTOOLS import SHMultiply
