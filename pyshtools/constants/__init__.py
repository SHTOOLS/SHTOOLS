"""
pyshtools constants.

This subpackage defines several constants used in analyzing gravity,
topography, and magnetic field data of the terrestrial planets. The constants
are organized by planet, and each object is an astropy Constant that possesses
the attributes name, value, unit, uncertainty, and reference. These constants
can be used in arithmetic operations with objects of the astropy class
Quantity.

Examples

Calculate the gravitational acceleration on the surface of Mars and
then to convert this to mGals:

    >>> Mars.gm / Mars.r**2
    <Quantity 3.7278663 m / s2>
    >>> (Mars.gm / Mars.r**2).to_value('mGal')
    372786.6303857397

Inspect a constant using the print function:

    >>> print(G)
      Name   = Gravitational constant
      Value  = 6.6743e-11
      Uncertainty  = 1.5e-15
      Unit  = m3 / (kg s2)
      Reference = CODATA 2018
"""

try:
    from astropy.constants import Constant
    from astropy.units.quantity import Quantity
except ImportError:
    raise ImportError('To use the pyshtools constant subpackage, you must '
                      'install astropy.')


# == Fundamental constants ==

from astropy.constants import G
from astropy.constants import mu0
from astropy.constants import codata

# == Constants organized by planet ===

from . import Mercury
from . import Venus
from . import Earth
from . import Moon
from . import Mars

# === Define __all__ ===

__all__ = ['Constant', 'Quantity', 'G', 'mu0', 'codata', 'Mercury', 'Venus',
           'Earth', 'Moon', 'Mars']
