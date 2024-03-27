"""
pyshtools constants.

This subpackage defines several constants used in analyzing gravity,
topography, and magnetic field data of the terrestrial planets. The constants
are organized by planet, and each object is an astropy Constant that possesses
the attributes name, value, unit, uncertainty, and reference. These constants
can be used in arithmetic operations with objects of the astropy class
Quantity.

Examples

Inspect a constant using the print function:

    >>> print(G)
      Name   = Gravitational constant
      Value  = 6.6743e-11
      Uncertainty  = 1.5e-15
      Unit  = m3 / (kg s2)
      Reference = CODATA 2018

Convert the orbital period of Callisto to days

    >>> Callisto.rotational_period.to('day')
    <Quantity 16.68901797 d>

Calculate the gravitational acceleration on the mean planetary radius of
Mercury and then return a value converted to mGals:

    >>> (Mercury.gm / Mercury.mean_radius**2).to_value('mGal')
    370218.70697392424
"""

try:
    from astropy.constants import codata
except ImportError:
    raise ImportError('To use the pyshtools constant subpackage, you must '
                      'install astropy.')


# == Fundamental constants ==

from astropy.constants import G
from astropy.constants import mu0
from astropy.constants import au

# == Constants organized by planet ===

from . import Sun
from . import Mercury
from . import Venus
from . import Earth
from . import Moon
from . import Mars
from . import Ceres
from . import Vesta
from . import Eros
from . import Jupiter
from . import Io
from . import Europa
from . import Ganymede
from . import Callisto
from . import Saturn
from . import Titan
from . import Enceladus
from . import Uranus
from . import Neptune
from . import Pluto
from . import Charon


# === Define __all__ ===

__all__ = ['G', 'mu0', 'codata', 'au', 'Sun', 'Mercury', 'Venus', 'Earth',
           'Moon', 'Mars', 'Vesta', 'Ceres', 'Eros', 'Jupiter', 'Io', 'Europa',
           'Ganymede', 'Callisto', 'Saturn', 'Titan', 'Enceladus', 'Uranus',
           'Neptune', 'Pluto', 'Charon']
