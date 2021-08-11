"""
pyshtools constants for the planet Mercury.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G


gm = _Constant(
    abbrev='gm_mercury',
    name='Gravitational constant times the mass of Mercury',
    value=2.2031815411154894e+13,
    unit='m3 / s2',
    uncertainty=1.9361909444154922e+5,
    reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
    'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
    '(2014), The gravity field, orientation, and ephemeris of Mercury '
    'from MESSENGER observations after three years in orbit, J. Geophys. '
    'Res. Planets, 119, 2417-2436, doi:10.1002/2014JE004675.')

mass = _Constant(
    abbrev='mass_mercury',
    name='Mass of Mercury',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_mercury and G.')

mean_radius = _Constant(
    abbrev='r_mercury',
    name='Mean radius of Mercury',
    value=2439.40197456433e3,
    unit='m',
    uncertainty=0.0,
    reference='gtmes_150v05: Smith, D. E., M. T. Zuber, R. J. Phillips, '
    'S. C. Solomon, S. A. Hauck II, F. G. Lemoine, E. Mazarico, G. A. '
    'Neumann, S. J. Peale, J.-L. Margot, C. L. Johnson, M. H. Torrence, '
    'M. E. Perry, D. D. Rowlands, S. Goossens, J. W. Head, A. H. Taylor '
    '(2012). Gravity field and internal structure of Mercury from '
    'MESSENGER. Science, 336, 214-217, doi:10.1126/science.1218809.')

r = mean_radius

density = _Constant(
    abbrev='density_mercury',
    name='Mean density of Mercury',
    value=3 * mass.value / (_np.pi * 4 * mean_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * mean_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         mean_radius.uncertainty /
                         (_np.pi * 4 * mean_radius.value**4))**2
                         ),
    reference='Derived from mass_mercury and r_mercury.')

g0 = _Constant(
    abbrev='g0_mercury',
    name='Surface gravity of Mercury, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_mercury and r_mercury.')

omega = _Constant(
    abbrev='omega_mercury',
    name='Angular spin rate of Mercury',
    value=6.1385108 * 2 * _np.pi / 360 / (24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.0,
    reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
    'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
    '(2014), The gravity field, orientation, and ephemeris of Mercury '
    'from MESSENGER observations after three years in orbit, J. Geophys. '
    'Res. Planets, 119, 2417-2436, doi:10.1002/2014JE004675.')

omega_orbit = _Constant(
    abbrev='omega_orbit_mercury',
    name='Angular rotation rate of Mercury about the Sun',
    value=2 * _np.pi / (87.969216879 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=6 * 2 * _np.pi / (87.969216879 * 24 * 60 * 60)**2,
    reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
    'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
    '(2014), The gravity field, orientation, and ephemeris of Mercury '
    'from MESSENGER observations after three years in orbit, J. Geophys. '
    'Res. Planets, 119, 2417-2436, doi:10.1002/2014JE004675.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'density', 'g0', 'omega',
           'omega_orbit']
