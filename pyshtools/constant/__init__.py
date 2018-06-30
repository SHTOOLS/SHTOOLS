"""
pyshtools constants.

This subpackage defines several constants used in analyzing gravity,
topography, and magnetic field data of the terrestrial planets. Each object is
an astropy Constant that possesses the attributes name, value, unit,
uncertainty, and reference. These constants can be used in arithmetic
operations with objects of the astropy class Quantity.

Examples

Calculate the gravitational acceleration on the surface of Mars and
then to convert this to mGals:

    >>> gm_mars / r_mars**2
    <Quantity 3.7278663 m / s2>
    >>> (gm_mars / r_mars**2).to_value('mGal')
    372786.6303857397

Inspect a constant using the print function:

    >>> print(G)
      Name   = Gravitational constant
      Value  = 6.67408e-11
      Uncertainty  = 3.1e-15
      Unit  = m3 / (kg s2)
      Reference = CODATA 2014
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np

try:
    from astropy.constants import Constant
    from astropy.units.quantity import Quantity
except ImportError:
    raise ImportError('To use the pyshtools constant subpackage, you must '
                      'install astropy.')


# == Fundamental constants ==

from astropy.constants import G
from astropy.constants import mu0

# == Mercury ===

from .Mercury import gm_mercury
from .Mercury import mass_mercury
from .Mercury import r_mercury
from .Mercury import density_mercury
from .Mercury import g0_mercury
from .Mercury import omega_mercury
from .Mercury import omega_orbit_mercury

# == Venus ==

from .Venus import gm_venus
from .Venus import mass_venus
from .Venus import r_venus
from .Venus import density_venus
from .Venus import g0_venus
from .Venus import omega_venus

# == Earth ==

from .Earth import gm_egm2008
from .Earth import mass_egm2008
from .Earth import omega_egm2008
from .Earth import a_wgs84
from .Earth import f_wgs84
from .Earth import gm_wgs84
from .Earth import mass_wgs84
from .Earth import omega_wgs84
from .Earth import gma_wgs84
from .Earth import b_wgs84
from .Earth import r3_wgs84
from .Earth import u0_wgs84

# == Moon ==

from .Moon import gm_moon
from .Moon import mass_moon
from .Moon import r_moon
from .Moon import density_moon
from .Moon import g0_moon
from .Moon import a_orbit_moon
from .Moon import omega_moon
from .Moon import i_solid_moon
from .Moon import beta_moon
from .Moon import gamma_moon

# == Mars ==

from .Mars import gm_mars
from .Mars import mass_mars
from .Mars import r_mars
from .Mars import density_mars
from .Mars import g0_mars
from .Mars import omega_mars
from .Mars import a_mars
from .Mars import b_mars
from .Mars import f_mars
from .Mars import u0_mars


# === Define groups of constants and __all__ ===

_constants_fundamental = ['G', 'mu0']

_constants_mercury = ['gm_mercury', 'mass_mercury', 'r_mercury',
                      'density_mercury', 'g0_mercury', 'omega_mercury',
                      'omega_orbit_mercury']

_constants_venus = ['gm_venus', 'mass_venus', 'r_venus', 'density_venus',
                    'g0_venus', 'omega_venus']

_constants_earth = ['gm_egm2008', 'mass_egm2008', 'omega_egm2008', 'a_wgs84',
                    'f_wgs84', 'gm_wgs84', 'mass_wgs84', 'omega_wgs84',
                    'gma_wgs84', 'b_wgs84', 'r3_wgs84', 'u0_wgs84']

_constants_moon = ['gm_moon', 'mass_moon', 'r_moon', 'density_moon',
                   'g0_moon', 'a_orbit_moon', 'omega_moon', 'i_solid_moon',
                   'beta_moon', 'gamma_moon']

_constants_mars = ['gm_mars', 'mass_mars', 'r_mars', 'density_mars', 'g0_mars',
                   'omega_mars', 'a_mars', 'b_mars', 'f_mars', 'u0_mars']

_constants = _constants_fundamental + _constants_mercury + _constants_venus \
             + _constants_earth + _constants_moon + _constants_mars

__all__ = ['Constant', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars'] \
           + _constants


# === Update doc string to list all constants with short descriptions ===

lines = ['\nThe following constants are available:\n',
         22 * '=' + ' ' + 58 * '=',
         'Name                   Description',
         22 * '=' + ' ' + 58 * '=']
sep = 22*'-' + ' ' + 58 * '-'

for module_name, module_list in [
        ('Fundamental constants', _constants_fundamental),
        ('Mercury', _constants_mercury),
        ('Venus', _constants_venus),
        ('Earth', _constants_earth),
        ('Moon', _constants_moon),
        ('Mars', _constants_mars)]:
    if module_name != 'Fundamental constants':
        lines.append(sep)
    lines.append('{0:22}'.format(module_name))
    lines.append(sep)
    for c in module_list:
        name = locals()[c].name
        if len(name) > 58:
            name = name[:55] + '...'
        lines.append('{0:22} {1}'.format(
            locals()[c].abbrev, name))

lines.append(lines[1])

__doc__ += '\n'.join(lines)


# == Clean up namespace ===

del lines
del sep
del c
del name
del module_name
del module_list
