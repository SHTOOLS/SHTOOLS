"""
pyshtools constants.

This subpackage defines several constants used in analyzing gravity,
topography, and magnetic field data of the terrestrial planets. Each object is
an astropy Constant that possesses the attributes name, value, unit,
uncertainty, and reference.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np

try:
    from astropy.constants import Constant

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

# === Update doc string to list all constants and short descriptions ===

_lines = ['\nThe following constants are available:\n',
          22 * '=' + ' ' + 58 * '=', 'Name                   Description',
          22 * '=' + ' ' + 58 * '=']
_sep = 22*'-' + ' ' + 58 * '-'

_lines.append('{0:22}'.format('Fundamental constants'))
_lines.append(_sep)
for _c in _constants_fundamental:
    _name = locals()[_c].name
    if len(_name) > 58:
        _name = _name[:55] + '...'
    _lines.append('{0:22} {1}'.format(
        locals()[_c].abbrev, _name))

_lines.append(_sep)
_lines.append('{0:22}'.format('Mercury'))
_lines.append(_sep)
for _c in _constants_mercury:
    _name = locals()[_c].name
    if len(_name) > 58:
        _name = _name[:55] + '...'
    _lines.append('{0:22} {1}'.format(
        locals()[_c].abbrev, _name))

_lines.append(_sep)
_lines.append('{0:22}'.format('Venus'))
_lines.append(_sep)
for _c in _constants_venus:
    _name = locals()[_c].name
    if len(_name) > 58:
        _name = _name[:55] + '...'
    _lines.append('{0:22} {1}'.format(
        locals()[_c].abbrev, _name))

_lines.append(_sep)
_lines.append('{0:22}'.format('Earth'))
_lines.append(_sep)
for _c in _constants_earth:
    _name = locals()[_c].name
    if len(_name) > 58:
        _name = _name[:55] + '...'
    _lines.append('{0:22} {1}'.format(
        locals()[_c].abbrev, _name))

_lines.append(_sep)
_lines.append('{0:22}'.format('The Moon'))
_lines.append(_sep)
for _c in _constants_moon:
    _name = locals()[_c].name
    if len(_name) > 58:
        _name = _name[:55] + '...'
    _lines.append('{0:22} {1}'.format(
        locals()[_c].abbrev, _name))

_lines.append(_sep)
_lines.append('{0:22}'.format('Mars'))
_lines.append(_sep)
for _c in _constants_mars:
    _name = locals()[_c].name
    if len(_name) > 58:
        _name = _name[:55] + '...'
    _lines.append('{0:22} {1}'.format(
        locals()[_c].abbrev, _name))

_lines.append(_lines[1])

__doc__ += '\n'.join(_lines)

del _lines
del _sep
