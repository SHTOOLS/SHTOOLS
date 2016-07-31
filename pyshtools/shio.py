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


# ---------------------------------------------------------------------
# --- Define a Python function that replaces the Fortran
# --- equivalent that uses different indexing conventions.
# ---------------------------------------------------------------------
def YilmIndexVector(i, l, m):
    return l**2 + (i - 1) * l + m
