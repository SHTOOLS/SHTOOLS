"""
pyshtools constants for the planet Uranus.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G
from astropy.constants import au as _au
from . import Sun as _Sun


gm = _Constant(
    abbrev='gm_uranus',
    name='Gravitational constant times the mass of the Uranus system',
    value=5794556.4e9,
    unit='m3 / s2',
    uncertainty=4.3e9,
    reference='Jacobson, R. A. (2014). The orbits of the Uranian satellites '
    'and rings, the gravity field of the Uranian system, and the orientation '
    'of the pole of Uranus. The Astronomical Journal, 148(5), 76. '
    'https://doi.org/10.1088/0004-6256/148/5/76')

mass = _Constant(
    abbrev='mass_uranus',
    name='Mass of Uranus',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_uranus and G.')

angular_velocity = _Constant(
    abbrev='angular_velocity_uranus',
    name='Angular spin rate of Uranus',
    value=2. * _np.pi / (17.239 * 60 * 60),
    unit='rad / s',
    uncertainty=0.009 * 2 * _np.pi / 17.239**2 / 60 / 60,
    reference='Desch, M. D., Connerney, J. E. P., & Kaiser, M. L. (1986). '
    'The rotation period of Uranus. Nature, 322(6074), 42–43. '
    'https://doi.org/10.1038/322042a0')

rotational_period = _Constant(
    abbrev='rotational_period_uranus',
    name='Rotational period of Uranus',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_uranus')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_uranus',
    name='Semimajor axis of the orbit of Uranus about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=19.18916464,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_uranus',
    name='Eccentricity of the orbit of Uranus about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=0.04725744,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_uranus',
    name='Inclination of the orbit of Uranus about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=0.77263783,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_uranus',
    name='Orbital angular velocity of Uranus',
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
    reference="Approximated using Kepler's third law, gm_sun, gm_uranus and "
    'orbit_semimajor_axis_uranus')

orbit_period = _Constant(
    abbrev='orbit_period_uranus',
    name='Orbital period of Uranus',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_uranus')

__all__ = ['gm', 'mass', 'angular_velocity', 'otational_period',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_period']
