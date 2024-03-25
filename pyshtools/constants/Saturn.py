"""
pyshtools constants for Saturn.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G
from astropy.constants import au as _au
from . import Sun as _Sun


orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_saturn',
    name='Semimajor axis of the orbit of Saturn about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=9.53667594,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_saturn',
    name='Eccentricity of the orbit of Saturn about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=0.05386179,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_saturn',
    name='Inclination of the orbit of Saturn about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=2.48599187,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

__all__ = ['orbit_semimajor_axis', 'orbit_eccentricity', 'orbit_inclination']
