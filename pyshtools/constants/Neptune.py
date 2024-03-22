"""
pyshtools constants for Neptune.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_neptune',
    name='Semimajor axis of the orbit of Neptune about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=30.06992276,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_neptune',
    name='Eccentricity of the orbit of Neptune about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=0.00859048,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_neptune',
    name='Inclination of the orbit of Neptune about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=1.77004347,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

__all__ = ['orbit_semimajor_axis', 'orbit_eccentricity', 'orbit_inclination']
