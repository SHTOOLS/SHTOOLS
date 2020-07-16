"""
Class interface for pyshtools.

pyshtools defines several classes that facilitate the interactive examination
of geographical gridded data and their associated spherical harmonic
coefficients. Superclasses are used to implement interface functions and
documentation whereas subclasses are used to handle different internal data
types.

Class structure:

    SHCoeffs
        SHRealCoeffs
        SHComplexCoeffs

    SHGrid
        DHRealGrid
        DHComplexGrid
        GLQRealGrid
        GLQComplexGrid

    SHWindow
        SHWindowCap
        SHWindowMask

    Slepian
        SlepianCap
        SlepianMask

    SlepianCoeffs

    SHGravCoeffs
        SHGravRealCoeffs

    SHMagCoeffs
        SHMagRealCoeffs

    SHGravGrid

    SHMagGrid

    SHTensor
        SHGravTensor
        SHMagTensor

    SHGeoid

For more information, see the documentation for the top level classes.
"""
from .shcoeffs import SHCoeffs
from .shcoeffs import SHRealCoeffs
from .shcoeffs import SHComplexCoeffs

from .shgrid import SHGrid
from .shgrid import DHRealGrid
from .shgrid import DHComplexGrid
from .shgrid import GLQRealGrid
from .shgrid import GLQComplexGrid

from .shwindow import SHWindow
from .shwindow import SHWindowCap
from .shwindow import SHWindowMask

from .slepian import Slepian
from .slepian import SlepianCap
from .slepian import SlepianMask

from .slepiancoeffs import SlepianCoeffs

from .shgravcoeffs import SHGravCoeffs
from .shgravcoeffs import SHGravRealCoeffs
from .shgravgrid import SHGravGrid
from .shtensor import SHGravTensor
from .shgeoid import SHGeoid

from .shmagcoeffs import SHMagCoeffs
from .shmagcoeffs import SHMagRealCoeffs
from .shmaggrid import SHMagGrid
from .shtensor import SHMagTensor


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['SHCoeffs', 'SHRealCoeffs', 'SHComplexCoeffs', 'SHGrid',
           'DHRealGrid', 'DHComplexGrid', 'GLQRealGrid', 'GLQComplexGrid',
           'SHWindow', 'SHWindowCap', 'SHWindowMask', 'Slepian', 'SlepianCap',
           'SlepianMask', 'SlepianCoeffs', 'SHGravCoeffs', 'SHGravRealCoeffs',
           'SHGravGrid', 'SHGravTensor', 'SHGeoid', 'SHMagCoeffs',
           'SHMagRealCoeffs', 'SHMagGrid', 'SHMagTensor']
