"""
pyshtools constants for the Sun.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_sun',
    name='Gravitational constant times the mass of the Sun',
    value=132712440041.279419e9,
    unit='m3 / s2',
    uncertainty=0.,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414.')

mass = _Constant(
    abbrev='mass_sun',
    name='Mass of the Sun',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_sun and G.')


__all__ = ['gm', 'mass']
