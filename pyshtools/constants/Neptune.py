"""
pyshtools constants for the planet Neptune.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G
from astropy.constants import au as _au
from . import Sun as _Sun


gm = _Constant(
    abbrev='gm_neptune',
    name='Gravitational constant times the mass of the Neptune system',
    value=6836527.100580397e9,
    unit='m3 / s2',
    uncertainty=10.e9,
    reference='Jacobson, R. A. (2009). The orbits of the Neptunian satellites '
    'and the orientation of the pole of Neptune. The Astronomical Journal, '
    '137(5), 4322. https://doi.org/10.1088/0004-6256/137/5/4322')

mass = _Constant(
    abbrev='mass_neptune',
    name='Mass of Neptune',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_neptune and G.')

angular_velocity = _Constant(
    abbrev='angular_velocity_neptune',
    name='Angular spin rate of Neptune',
    value=536.3128492 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=1.66 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    reference='Warwick, J. W., et al. (1989). Voyager Planetary Radio '
    'Astronomy at Neptune. Science, 246(4936), 1498–1501. '
    'https://doi.org/10.1126/science.246.4936.1498')

rotational_period = _Constant(
    abbrev='rotational_period_neptune',
    name='Rotational period of Neptune',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_neptune')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_neptune',
    name='Semimajor axis of the orbit of Neptune about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
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
    'the mean ecliptic and equinox of J2000',
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
    'the mean ecliptic and equinox of J2000',
    value=1.77004347,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_neptune',
    name='Orbital angular velocity of Neptune',
    value=_np.sqrt((_Sun.gm.value + gm.value) /
                   (_au.value * orbit_semimajor_axis.value)**3),
    unit='rad / s',
    uncertainty=_np.sqrt(
        _Sun.gm.uncertainty**2 / 4. / (_Sun.gm.value + gm.value) /
        (_au.value * orbit_semimajor_axis.value)**3 +
        gm.uncertainty**2 / 4. / (_Sun.gm.value + gm.value) /
        (_au.value * orbit_semimajor_axis.value)**3 +
        9. * (_au.value * orbit_semimajor_axis.uncertainty)**2 *
        (_Sun.gm.value + gm.value) / 4. /
        (_au.value * orbit_semimajor_axis.value)**5),
    reference="Approximated using Kepler's third law, gm_sun, gm_neptune and "
    'orbit_semimajor_axis_neptune')

orbit_period = _Constant(
    abbrev='orbit_period_neptune',
    name='Orbital period of Neptune',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_neptune')

__all__ = ['gm', 'mass', 'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_period']
