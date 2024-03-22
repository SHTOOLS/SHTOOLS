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
    abbrev='mean_radius_mercury',
    name='Mean radius of Mercury',
    value=2439472.7,
    unit='m',
    uncertainty=0.0,
    reference='USGS_SPG_shape: Maia, J. (2024). Spherical harmonic models of '
    'the shape of Mercury [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.10809345')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_mercury',
    name='Volume equivalent radius of Mercury',
    value=2439473.1,
    unit='m',
    uncertainty=0.,
    reference='Computed using USGS_SPG_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_mercury',
    name='Volume of Mercury',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_mercury')

mean_density = _Constant(
    abbrev='mean_density_mercury',
    name='Mean density of Mercury',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_mercury and '
    'volume_equivalent_radius_mercury.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_mercury',
    name='Gravity at the mean radius of Mercury, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_mercury and mean_radius_mercury.')

angular_velocity = _Constant(
    abbrev='angular_velocity_mercury',
    name='Angular spin rate of Mercury',
    value=6.1385108 * 2 * _np.pi / 360 / (24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.0,
    reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
    'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
    '(2014), The gravity field, orientation, and ephemeris of Mercury '
    'from MESSENGER observations after three years in orbit, J. Geophys. '
    'Res. Planets, 119, 2417-2436, doi:10.1002/2014JE004675.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_mercury',
    name='Angular rotation rate of Mercury about the Sun',
    value=2 * _np.pi / (87.969216879 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=6 * 2 * _np.pi / (87.969216879 * 24 * 60 * 60)**2,
    reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
    'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
    '(2014), The gravity field, orientation, and ephemeris of Mercury '
    'from MESSENGER observations after three years in orbit, J. Geophys. '
    'Res. Planets, 119, 2417-2436, doi:10.1002/2014JE004675.')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_mercury',
    name='Semimajor axis of the orbit of Mercury about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=0.38709927,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_mercury',
    name='Eccentricity of the orbit of Mercury about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=0.20563593,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_mercury',
    name='Inclination of the orbit of Mercury about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=7.00497902,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'mean_density', 'gravity_mean_radius', 'angular_velocity',
           'orbit_angular_velocity', 'orbit_semimajor_axis',
           'orbit_eccentricity', 'orbit_inclination']
