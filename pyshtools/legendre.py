"""
pyshtools Legendre Functions.

This submodule of pyshtools defines the following functions:

"Geodesy" 4-pi normalized
-----------------------
PlmBar        Compute all the geodesy-normalized associated Legendre functions.
PlmBar_d1     Compute all the geodesy-normalized associated Legendre functions
              and first derivatives.
PlBar         Compute all the geodesy-normalized Legendre polynomials.
PlBar_d1      Compute all the geodesy-normalized Legendre Polynomials and first
              derivatives.

Orthonormalized
---------------
PlmON         Compute all the orthonormalized associated Legendre functions.
PlmON_d1      Compute all the orthonormalized associated Legendre functions and
              first derivatives.
PlON          Compute all the orthonormalized Legendre polynomials.
PlON_d1       Compute all the orthonormalized Legendre polynomials and first
              derivatives.

Schmidt normalized
------------------
PlmSchmidt    Compute all the Schmidt-normalized associated Legendre functions.
PlmSchmidt_d1 Compute all the Schmidt-normalized associated Legendre functions
              and first derivatives.
PlSchmidt     Compute all the Schmidt-normalized Legendre polynomials.
PlSchmidt_d1  Compute all the Schmidt-normalized Legendre polynomials and first
              derivatives.

Unnormalized
------------
PLegendreA    Compute all the unnormalized associated Legendre functions.
PLegendreA_d1 Compute all the unnormalized associated Legendre functions and
              first derivatives.
PLegendre     Compute all the unnormalized Legendre polynomials.
PLegendre_d1  Compute all the unnormalized Legendre polynomials and first
              derivatives.

Other
-----
PlmIndex      Compute the index of an array of Legendre function corresponding
              to degree L and angular order M.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

from ._SHTOOLS import PlmBar
from ._SHTOOLS import PlmBar_d1
from ._SHTOOLS import PlBar
from ._SHTOOLS import PlBar_d1
from ._SHTOOLS import PlmON
from ._SHTOOLS import PlmON_d1
from ._SHTOOLS import PlON
from ._SHTOOLS import PlON_d1
from ._SHTOOLS import PlmSchmidt
from ._SHTOOLS import PlmSchmidt_d1
from ._SHTOOLS import PlSchmidt
from ._SHTOOLS import PlSchmidt_d1
from ._SHTOOLS import PLegendreA
from ._SHTOOLS import PLegendreA_d1
from ._SHTOOLS import PLegendre
from ._SHTOOLS import PLegendre_d1


# ---------------------------------------------------------------------
# --- Define a Python function that replaces the Fortran
# --- equivalent that uses different indexing conventions.
# ---------------------------------------------------------------------
def PlmIndex(l, m):
    return (l * (l + 1)) // 2 + m
