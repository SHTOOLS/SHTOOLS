"""
pyshtools utils.

This submodule of pyshtools defines the following functions:

MakeCircleCoord    Compute coordinates of a circle placed at a given latitude
                   and longitude.
MakeEllipseCoord   Compute coordinates of an ellipse placed at a given latitude
                   and longitude.
Wigner3j           Compute the Wigner-3j symbols for all allowable values of J.
DHaj               Compute the latitudinal weights used in the Driscoll and
                   Healy (1994) spherical harmonic transforms.
figstyle           Set Matplotlib parameters for creating publication quality
                   graphics.
"""
from ..backends.shtools import MakeCircleCoord
from ..backends.shtools import MakeEllipseCoord
from ..backends.shtools import Wigner3j
from ..backends.shtools import DHaj
from .figstyle import figstyle


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['MakeCircleCoord', 'MakeEllipseCoord', 'Wigner3j', 'DHaj',
           'figstyle']
