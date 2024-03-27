"""
pyshtools constants for the planet Earth.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
from types import SimpleNamespace as _SimpleNamespace
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G
from astropy.constants import au as _au
from . import Sun as _Sun
from . import Moon as _Moon


_a_wgs84 = _Constant(
    abbrev='a_wgs84',
    name='Semimajor axis of the WGS84 ellipsoid',
    value=6378137.0,
    unit='m',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_f_wgs84 = _Constant(
    abbrev='f_wgs84',
    name='Flattening of the WGS84 ellipsoid',
    value=1/298.257223563,
    unit='',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_gm_wgs84 = _Constant(
    abbrev='gm_wgs84',
    name='Gravitational constant times the mass of Earth for the '
         'WGS84 geodetic reference system, including the atmosphere',
    value=3986004.418e8,
    unit='m3 / s2',
    uncertainty=0.008e8,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_mass_wgs84 = _Constant(
    abbrev='mass_wgs84',
    name='Mass of Earth for the WGS84 geodetic reference system, '
    'including the atmosphere',
    value=_gm_wgs84.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((_gm_wgs84.uncertainty / _G.value)**2 +
                         (_gm_wgs84.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_wgs84 and G.')

_omega_wgs84 = _Constant(
    abbrev='omega_wgs84',
    name='Angular rotation rate of Earth for the WGS84 geodetic reference '
    'system',
    value=7292115.0e-11,
    unit='rad / s',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_gma_wgs84 = _Constant(
    abbrev='gma_wgs84',
    name="Gravitational constant times the mass of Earth's atmosphere "
    'for the WGS84 geodetic reference system',
    value=3.5e8,
    unit='m3 / s2',
    uncertainty=0.1e8,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_b_wgs84 = _Constant(
    abbrev='b_wgs84',
    name='Semiminor axis of the WGS84 ellipsoid',
    value=6356752.3142,
    unit='m',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_r3_wgs84 = _Constant(
    abbrev='r3_wgs84',
    name="Radius of a sphere with volume equal to that of the WGS84 "
    'ellipsoid',
    value=6371000.7900,
    unit='m',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_u0_wgs84 = _Constant(
    abbrev='u0_wgs84',
    name="Theoretical normal gravity potential of the WGS84 ellipsoid",
    value=62636851.7146,
    unit='m2 / s2',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

_rotational_period_wgs84 = _Constant(
    abbrev='rotational_period_wgs84',
    name='Rotational period of the WGS84 ellipsoid',
    value=2. * _np.pi / _omega_wgs84.value,
    unit='s',
    uncertainty=0.,
    reference='Derived from wgs84.omega')

wgs84 = _SimpleNamespace(__doc__='WGS84 geodetic reference system')
wgs84.a = _a_wgs84
wgs84.f = _f_wgs84
wgs84.gm = _gm_wgs84
wgs84.mass = _mass_wgs84
wgs84.omega = _omega_wgs84
wgs84.gma = _gma_wgs84
wgs84.b = _b_wgs84
wgs84.r3 = _r3_wgs84
wgs84.u0 = _u0_wgs84
wgs84.rotational_period = _rotational_period_wgs84

_gm_egm2008 = _Constant(
    abbrev='gm_egm2008',
    name='Gravitational constant times the mass of Earth for the model '
         'EGM2008, including the atmosphere',
    value=3986004.415e+8,
    unit='m3 / s2',
    uncertainty=0.0,
    reference='Pavlis N. K., S. A. Holmes, S. C. Kenyon, and J. K. Factor '
    '(2012). The development and evaluation of the Earth Gravitational '
    'Model 2008 (EGM2008). J. Geophys. Res., 117, B04406, '
    'doi:10.1029/2011JB008916.')

_mass_egm2008 = _Constant(
    abbrev='mass_egm2008',
    name='Mass of Earth for the model EGM2008, including the atmosphere',
    value=_gm_egm2008.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((_gm_egm2008.uncertainty / _G.value)**2 +
                         (_gm_egm2008.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_egm2008 and G.')

_omega_egm2008 = _Constant(
    abbrev='omega_egm2008',
    name='Angular spin rate of Earth for the model EGM2008',
    value=7292115.0e-11,
    unit='rad / s',
    uncertainty=0.0,
    reference='Pavlis N. K., S. A. Holmes, S. C. Kenyon, and J. K. Factor '
    '(2012). The development and evaluation of the Earth Gravitational '
    'Model 2008 (EGM2008). J. Geophys. Res., 117, B04406, '
    'doi:10.1029/2011JB008916.')

_rotational_period_egm2008 = _Constant(
    abbrev='rotational_period_egm2008',
    name='Rotational period of the egm2008 ellipsoid',
    value=2. * _np.pi / _omega_egm2008.value,
    unit='s',
    uncertainty=0.,
    reference='Derived from egm2008.omega')

egm2008 = _SimpleNamespace(__doc__='EGM2008')
egm2008.gm = _gm_egm2008
egm2008.mass = _mass_egm2008
egm2008.omega = _omega_egm2008
egm2008.rotational_period = _rotational_period_egm2008

dynamical_flattening_earth = _Constant(
    abbrev='dynamical_flattening_earth',
    name="Dynamical flattening of the Earth",
    value=3273795e-9,
    unit='',
    uncertainty=1e-9,
    reference='Petit, G., & Luzum, B. (2010). IERS conventions (2010) '
    '(No. IERS-TN-36). BUREAU INTERNATIONAL DES POIDS '
    'ET MESURES SEVRES (FRANCE).')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_earth_moon_barycenter',
    name='Semimajor axis of the orbit of the Earth-Moon barycenter about the '
    'Sun, with respect to the mean ecliptic and equinox of J2000',
    value=1.00000261,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_earth_moon_barycenter',
    name='Eccentricity of the orbit of the Earth-Moon barycenter about the '
    'Sun, with respect to the mean ecliptic and equinox of J2000',
    value=0.01671123,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_earth_moon_barycenter',
    name='Inclination of the orbit of the Earth-Moon barycenter about the '
    'Sun, with respect to the mean ecliptic and equinox of J2000',
    value=-0.00001531,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_earth_moon_barycenter',
    name='Orbital angular velocity of Earth-Moon barycenter',
    value=_np.sqrt((_Sun.gm.value + egm2008.gm.value + _Moon.gm.value) /
                   (_au.value * orbit_semimajor_axis.value)**3),
    unit='rad / s',
    uncertainty=_np.sqrt(
        _Sun.gm.uncertainty**2 / 4. /
        (_Sun.gm.value + egm2008.gm.value + _Moon.gm.value) /
        (_au.value * orbit_semimajor_axis.value)**3 +
        egm2008.gm.uncertainty**2 / 4. /
        (_Sun.gm.value + egm2008.gm.value + _Moon.gm.value) /
        (_au.value * orbit_semimajor_axis.value)**3 +
        _Moon.gm.uncertainty**2 / 4. /
        (_Sun.gm.value + egm2008.gm.value + _Moon.gm.value) /
        (_au.value * orbit_semimajor_axis.value)**3 +
        9. * (_au.value * orbit_semimajor_axis.uncertainty)**2 *
        (_Sun.gm.value + egm2008.gm.value + _Moon.gm.value) / 4. /
        (_au.value * orbit_semimajor_axis.value)**5),
    reference="Approximated using Kepler's third law, gm_sun, gm_earth, "
    'gm_moon and orbit_semimajor_axis_earth_moon_barycenter')

orbit_period = _Constant(
    abbrev='orbit_period_earth_moon_barycenter',
    name='Orbital period of the Earth-Moon barycenter',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_earth_moon_barycenter')

__all__ = ['wgs84', 'egm2008', 'dynamical_flattening_earth',
           'orbit_semimajor_axis', 'orbit_eccentricity', 'orbit_inclination',
           'orbit_angular_velocity', 'orbit_period']
