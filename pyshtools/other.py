"""
pyshtools Other Routines.

This submodule of pyshtools defines the following functions:

MakeCircleCoord    Compute coordinates of a circle placed at a given latitude
                   and longitude.
MakeEllipseCoord   Compute coordinates of an ellipse placed at a given latitude
                   and longitude.
RandomN            Return a pseudo uniform random deviate between 0 and 1 using
                   the algorithm of Park and Miller with a Marsaglia shift
                   sequence.
RandomGaussian     Return a pseudo Gaussian deviate of zero mean and unit
                   variance.
PreGLQ             Calculate the weights and nodes used in integrating a
                   function by Gauss-Legendre quadrature.
EigValVecSym       Compute the eigenvalues and eigenvectors of a real
                   symmetric matrix.
EigValVecSymTri    Compute the eigenvalues and eigenvectors of a real
                   symmetric tridiagonal matrix.
EigValSym          Compute the eigenvalues of a real symmetric matrix.
Wigner3j           Compute the Wigner-3j symbols for all allowable values of J.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

from ._SHTOOLS import MakeCircleCoord
from ._SHTOOLS import MakeEllipseCoord
from ._SHTOOLS import RandomN
from ._SHTOOLS import RandomGaussian
from ._SHTOOLS import PreGLQ
from ._SHTOOLS import EigValVecSym
from ._SHTOOLS import EigValVecSymTri
from ._SHTOOLS import EigValSym
from ._SHTOOLS import Wigner3j
