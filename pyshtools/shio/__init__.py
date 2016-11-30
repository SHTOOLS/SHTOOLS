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
read_icgem_gfc   Read spherical harmonic coefficients or associated errors
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

from ..shtools import SHRead
from ..shtools import SHReadH
from ..shtools import SHReadError
from ..shtools import SHReadErrorH
from ..shtools import SHRead2
from ..shtools import SHRead2Error
from ..shtools import SHReadJPL
from ..shtools import SHReadJPLError
from ..shtools import SHCilmToCindex
from ..shtools import SHCindexToCilm
from ..shtools import SHCilmToVector
from ..shtools import SHVectorToCilm
from ..shtools import SHrtoc
from ..shtools import SHctor

from .icgem import read_icgem_gfc
from .yilm_index_vector import YilmIndexVector
