"""
pyshtools Gravity and Magnetics Routines.

This subpackage of pyshtools defines the following functions:

Gravity routines
----------------
MakeGravGridDH      Create 2D cylindrical maps on a flattened and rotating
                    ellipsoid of all three components of the gravity field,
                    the gravity disturbance, and the gravitational potential.
MakeGravGridPoint   Compute the gravity disturbance at a single point.
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
MakeMagGridPoint    Compute the magnetic field at a single point.
MakeMagGradGridDH   Calculate the components of the magnetic field tensor
                    on a flattened ellipsoid.
mag_spectrum        Compute the spectrum of either the magnetic potential
                    or magnetic field strength.
"""
from ..backends.shtools import MakeGravGridDH
from ..backends.shtools import MakeGravGridPoint
from ..backends.shtools import MakeGravGradGridDH
from ..backends.shtools import MakeGeoidGridDH
from ..backends.shtools import CilmPlusDH
from ..backends.shtools import CilmMinusDH
from ..backends.shtools import CilmPlusRhoHDH
from ..backends.shtools import CilmMinusRhoHDH
from ..backends.shtools import BAtoHilmDH
from ..backends.shtools import BAtoHilmRhoHDH
from ..backends.shtools import DownContFilterMA
from ..backends.shtools import DownContFilterMC
from ..backends.shtools import NormalGravity
from ..backends.shtools import MakeMagGridDH
from ..backends.shtools import MakeMagGridPoint
from ..backends.shtools import MakeMagGradGridDH

from .mag_spectrum import mag_spectrum


__all__ = ['MakeGravGridDH', 'MakeGravGradGridDH', 'MakeGeoidGridDH',
           'CilmPlusDH', 'CilmMinusDH', 'CilmPlusRhoHDH', 'CilmMinusRhoHDH',
           'BAtoHilmDH', 'BAtoHilmRhoHDH', 'DownContFilterMA',
           'DownContFilterMC', 'NormalGravity', 'MakeMagGridDH',
           'MakeMagGradGridDH', 'mag_spectrum', 'MakeGravGridPoint',
           'MakeMagGridPoint']
