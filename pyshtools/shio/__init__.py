"""
pyshtools Spherical Harmonic I/O, Storage, and Conversion Routines.

This submodule of pyshtools defines the following functions:

Spherical harmonic I/O
----------------------
shread           Read shtools-formatted spherical harmonic coefficients from a
                 text file.
shwrite          Write shtools-formatted spherical harmonic coefficients to a
                 text file.
read_dov         Read spherical harmonic coefficients from a text file
                 formatted as [degree, order, value].
write_dov        Write spherical harmonic coefficients to a text file formatted
                 as [degree, order, value].
read_bshc        Read real spherical harmonic coefficients from a binary
                 bshc-formatted file.
write_bshc       Write real spherical harmonic coefficients to a binary
                 bshc-formatted file.
read_icgem_gfc   Read real spherical harmonic gravitational potential
                 coefficients and associated errors from an ICGEM GFC formatted
                 file.
write_icgem_gfc  Write real spherical harmonic gravitational potential
                 coefficients and associated errors to an ICGEM GFC formatted
                 file.
read_igrf        Read IGRF real spherical harmonic coefficients, and return the
                 magnetic potential coefficients for the specified year.
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
convert          Convert an array of spherical harmonic coefficients to a
                 different normalization.
SHrtoc           Convert real spherical harmonics to complex form.
SHctor           Convert complex spherical harmonics to real form.
"""
from ..backends.shtools import SHRead2
from ..backends.shtools import SHRead2Error
from ..backends.shtools import SHReadJPL
from ..backends.shtools import SHReadJPLError
from ..backends.shtools import SHCilmToCindex
from ..backends.shtools import SHCindexToCilm
from ..backends.shtools import SHCilmToVector
from ..backends.shtools import SHVectorToCilm
from ..backends.shtools import SHrtoc
from ..backends.shtools import SHctor

from .convert import convert
from .shtools import shread
from .shtools import shwrite
from .dov import read_dov
from .dov import write_dov
from .icgem import read_icgem_gfc
from .icgem import write_icgem_gfc
from .read_igrf import read_igrf
from .bshc import read_bshc
from .bshc import write_bshc
from .yilm_index_vector import YilmIndexVector

del shtools  # noqa: F821
del dov  # noqa: F821
del bshc  # noqa: F821
del icgem  # noqa: F821
del yilm_index_vector  # noqa: F821


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['SHRead2', 'SHRead2Error', 'SHReadJPL', 'SHReadJPLError',
           'SHCilmToCindex', 'SHCindexToCilm', 'SHCilmToVector',
           'SHVectorToCilm', 'SHrtoc', 'SHctor', 'convert', 'shread',
           'shwrite', 'read_dov', 'write_dov', 'read_bshc', 'write_bshc',
           'read_icgem_gfc', 'write_icgem_gfc', 'read_igrf', 'YilmIndexVector']
