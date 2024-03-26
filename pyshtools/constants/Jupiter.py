"""
pyshtools constants for the planet Jupiter.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G
from astropy.constants import au as _au
from . import Sun as _Sun


gm = _Constant(
    abbrev='gm_jupiter',
    name='Gravitational constant times the mass of Jupiter',
    value=1.266432336e+17,
    unit='m3 / s2',
    uncertainty=0.,
    reference='Kaspi, Y., Galanti, E., Park, R. S., Duer, K., Gavriel, N., '
    'Durante, D., Iess, L., Parisi, M., Buccino, D. R., Guillot, T., '
    'Stevenson, D. J., & Bolton, S. J. (2023). Observational evidence for '
    'cylindrically oriented zonal flows on Jupiter. Nature Astronomy, 7(12), '
    '1463–1472. https://doi.org/10.1038/s41550-023-02077-8')

mass = _Constant(
    abbrev='mass_jupiter',
    name='Mass of Jupiter',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_jupiter and G.')

angular_velocity = _Constant(
    abbrev='angular_velocity_jupiter',
    name='Angular spin rate of Jupiter (system III)',
    value=0.00017585327063385654,
    unit='rad / s',
    uncertainty=1.4765300375888075e-11,
    reference='Yu, Z. J., & Russell, C. T. (2009). Rotation period of Jupiter '
    'from the observation of its magnetic field. Geophysical Research '
    'Letters, 36(20). https://doi.org/10.1029/2009GL040094')

rotational_period = _Constant(
    abbrev='rotational_period_jupiter',
    name='Rotational period of Jupiter (system III)',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_jupiter')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_jupiter',
    name='Semimajor axis of the orbit of Jupiter about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=5.20288700,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_jupiter',
    name='Eccentricity of the orbit of Jupiter about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=0.04838624,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_jupiter',
    name='Inclination of the orbit of Jupiter about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=1.30439695,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_jupiter',
    name='Orbital angular velocity of Jupiter',
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
    reference="Approximated using Kepler's third law, gm_sun, gm_jupiter and "
    'orbit_semimajor_axis_jupiter')

orbit_period = _Constant(
    abbrev='orbit_period_jupiter',
    name='Orbital period of Jupiter',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_uranus')

__all__ = ['gm', 'mass', 'angular_velocity', 'rotational_period',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_period']
