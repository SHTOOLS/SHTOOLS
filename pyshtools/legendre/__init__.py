"""
pyshtools Legendre Functions.

This subpackage of pyshtools defines the following functions:

Convenience functions
---------------------
legendre      Compute all the associated Legendre functions up to a maximum
              degree and order.
legendre_lm   Compute the associated Legendre function for a specific degree l
              and order m.

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
from ..backends.shtools import PlmBar
from ..backends.shtools import PlmBar_d1
from ..backends.shtools import PlBar
from ..backends.shtools import PlBar_d1
from ..backends.shtools import PlmON
from ..backends.shtools import PlmON_d1
from ..backends.shtools import PlON
from ..backends.shtools import PlON_d1
from ..backends.shtools import PlmSchmidt
from ..backends.shtools import PlmSchmidt_d1
from ..backends.shtools import PlSchmidt
from ..backends.shtools import PlSchmidt_d1
from ..backends.shtools import PLegendreA
from ..backends.shtools import PLegendreA_d1
from ..backends.shtools import PLegendre
from ..backends.shtools import PLegendre_d1

from .legendre_functions import legendre
from .legendre_functions import legendre_lm
from .plm_index import PlmIndex

del legendre_functions  # noqa: F821


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['PlmBar', 'PlmBar_d1', 'PlBar', 'PlBar_d1', 'PlmON', 'PlmON_d1',
           'PlON', 'PlON_d1', 'PlmSchmidt', 'PlmSchmidt_d1', 'PlSchmidt',
           'PlSchmidt_d1', 'PLegendreA', 'PLegendreA_d1', 'PLegendre',
           'PLegendre_d1', 'legendre', 'legendre_lm', 'PlmIndex']
