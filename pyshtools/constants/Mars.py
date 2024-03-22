"""
pyshtools constants for the planet Mars.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_mars',
    name='Gravitational constant times the mass of Mars',
    value=0.4282837581575610e+14,
    unit='m3 / s2',
    uncertainty=0.18167460e+6,
    reference='Konopliv A. S., R. S. Park, W. M. Folkner (2016). '
    'An improved JPL Mars gravity field and orientation from Mars orbiter '
    'and lander tracking data, Icarus, 274, 253-260, '
    'doi:10.1016/j.icarus.2016.02.052')

mass = _Constant(
    abbrev='mass_mars',
    name='Mass of Mars',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_mars and G.')

mean_radius = _Constant(
    abbrev='mean_radius_mars',
    name='Mean radius of Mars',
    value=3389.5e3,
    unit='m',
    uncertainty=0.0,
    reference='MOLA_shape: Wieczorek, M. (2024). Spherical harmonic models of '
    'the shape of Mars (1.0.0) [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.10794059')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_mars',
    name='Volume equivalent radius of Mars',
    value=3389513.3,
    unit='m',
    uncertainty=0.,
    reference='Computed using MOLA_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_mars',
    name='Volume of Mars',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_mars')

mean_density = _Constant(
    abbrev='mean_density_mars',
    name='Mean density of Mars',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_mars and volume_equivalent_radius_mars.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_mars',
    name='Gravity at the mean radius of Mars, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_mars and mean_radius_mars.')

angular_velocity = _Constant(
    abbrev='angular_velocity_mars',
    name='Angular spin rate of Mars',
    value=350.891985307 * 2 * _np.pi / 360 / (24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.000000003 * 2 * _np.pi / 360 / (24 * 60 * 60),
    reference='Konopliv A. S., R. S. Park, W. M. Folkner (2016). '
    'An improved JPL Mars gravity field and orientation from Mars orbiter '
    'and lander tracking data, Icarus, 274, 253-260, '
    'doi:10.1016/j.icarus.2016.02.052')

a = _Constant(
    abbrev='semimajor_axis_mars',
    name='Semimajor axis of the Mars reference ellipsoid',
    value=3395428.0,
    unit='m',
    uncertainty=19.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

b = _Constant(
    abbrev='semiminor_axis_mars',
    name='Semiminor axis of the Mars reference ellipsoid',
    value=3377678.0,
    unit='m',
    uncertainty=19.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

f = _Constant(
    abbrev='flattening_mars',
    name='Flattening of the Mars reference ellipsoid',
    value=(a.value - b.value) / a.value,
    unit='',
    uncertainty=_np.sqrt((a.uncertainty * (a.value - b.value)
                         / a.value**2)**2
                         + (_np.sqrt(a.uncertainty**2 +
                            b.uncertainty**2) / a.value)**2
                         ),
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

u0 = _Constant(
    abbrev='normal_gravity_potential_mars',
    name='Theoretical normal gravity potential of the reference ellipsoid',
    value=12654875.0,
    unit='m2 / s2',
    uncertainty=69.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_mars',
    name='Semimajor axis of the orbit of Mars about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=1.52371034,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_mars',
    name='Eccentricity of the orbit of Mars about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=0.09339410,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_mars',
    name='Inclination of the orbit of Mars about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000, valid for the time interval '
    '1800-2050 AD',
    value=1.84969142,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'mean_density',
           'volume_equivalent_radius', 'volume', 'gravity_mean_radius',
           'angular_velocity', 'a', 'b', 'f', 'u0', 'orbit_semimajor_axis',
           'orbit_eccentricity', 'orbit_inclination']
