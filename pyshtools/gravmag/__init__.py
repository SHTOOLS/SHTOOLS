"""
pyshtools Gravity and Magnetics Routines.

This subpackage of pyshtools defines the following functions:

Gravity routines
----------------
MakeGravGridDH      Create 2D cylindrical maps on a flattened and rotating
                    ellipsoid of all three components of the gravity field,
                    the gravity disturbance, and the gravitational potential.
MakeGravGradGridDH  Calculate the components of the gravity "gradient" tensor
                    on a flattened ellipsoid.
MakeGeoidGridDH     Create a global map of the geoid.
CilmPlusDH          Calculate the gravitational potential exterior to relief
                    along a spherical interface using the finite-amplitude
                    algorithm of Wieczorek and Phillips (1998) on a Driscoll
                    and Healy (1994) grid.
CilmMinusDH         Calculate the gravitational potential interior to relief
                    along to a spherical interface using the finite-amplitude
                    algorithm of Wieczorek and Phillips (1998) on a Driscoll
                    and Healy (1994) grid.
CilmPlusRhoHDH      Calculate the gravitational potential exterior to relief
                    along a spherical interface with laterally varying density
                    using the finite amplitude algorithm of Wieczorek (2007)
                    on a Driscoll and Healy (1994) grid.
CilmMinusRhoHDH     Calculate the gravitational potential interior to relief
                    along a spherical interface with laterally varying density
                    using the finite amplitude algorithm of Wieczorek (2007)
                    on a Driscoll and Healy (1994) grid.
BAtoHilmDH          Calculate iteratively the relief along an interface with
                    constant density contrast that corresponds to a given
                    Bouguer anomaly using the algorithm of Wieczorek and
                    Phillips (1998).
BAtoHilmRhoHDH      Iteratively calculate the relief along an interface with
                    laterally varying density contrast that corresponds to a
                    given Bouguer anomaly using the algorithm of Wieczorek and
                    Phillips (1998).
DownContFilterMA    Compute the minimum-amplitude downward continuation filter
                    of Wieczorek and Phillips (1998).
DownContFilterMC    Calculate a minimum-curvature downward continuation filter
                    for a given spherical harmonic degree.
NormalGravity       Calculate the normal gravity on a flattened ellipsoid
                    using the formula of Somigliana.

Magnetics routines
------------------
MakeMagGridDH       Create 2D cylindrical maps on a flattened ellipsoid of
                    all three vector components of the magnetic field, the
                    magnitude of the magnetic field, and the magnetic
                    potential.
MakeMagGradGridDH   Calculate the components of the magnetic field tensor
                    on a flattened ellipsoid.
mag_spectrum        Compute the spectrum of either the magnetic potential
                    or magnetic field strength.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

from ..shtools import MakeGravGridDH
from ..shtools import MakeGravGradGridDH
from ..shtools import MakeGeoidGridDH
from ..shtools import CilmPlusDH
from ..shtools import CilmMinusDH
from ..shtools import CilmPlusRhoHDH
from ..shtools import CilmMinusRhoHDH
from ..shtools import BAtoHilmDH
from ..shtools import BAtoHilmRhoHDH
from ..shtools import DownContFilterMA
from ..shtools import DownContFilterMC
from ..shtools import NormalGravity
from ..shtools import MakeMagGridDH
from ..shtools import MakeMagGradGridDH

from .mag_spectrum import mag_spectrum
