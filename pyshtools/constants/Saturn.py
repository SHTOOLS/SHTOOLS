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


gm = _Constant(
    abbrev='gm_saturn',
    name='Gravitational constant times the mass of Saturn',
    value=37931206.234e9,
    unit='m3 / s2',
    uncertainty=0.726e9,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

mass = _Constant(
    abbrev='mass_saturn',
    name='Mass of Saturn',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_saturn and G.')

angular_velocity = _Constant(
    abbrev='angular_velocity_saturn',
    name='Angular spin rate of Saturn',
    value=818.1387763691 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=2.05 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    reference='Mankovich, C., Marley, M. S., Fortney, J. J., & Movshovitz, N. '
    "(2019). Cassini Ring Seismology as a Probe of Saturn's Interior. I. '
    'Rigid Rotation. The Astrophysical Journal, 871(1), 1. '
    'https://doi.org/10.3847/1538-4357/aaf798')

rotational_period = _Constant(
    abbrev='rotational_period_saturn',
    name='Rotational period of Saturn',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_saturn')

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

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_saturn',
    name='Orbital angular velocity of saturn',
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
    reference="Approximated using Kepler's third law, gm_sun, gm_saturn and "
    'orbit_semimajor_axis_saturn')

orbit_period = _Constant(
    abbrev='orbit_period_saturn',
    name='Orbital period of Saturn',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_saturn')

__all__ = ['gm', 'mass', 'angular_velocity', 'rotational_period',
           'orbit_semimajor_axis', 'orbit_eccentricity', 'orbit_inclination',
           'orbit_angular_velocity', 'orbit_period']
