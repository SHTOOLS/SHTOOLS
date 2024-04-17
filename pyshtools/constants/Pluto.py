"""
pyshtools constants for Pluto.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G
from astropy.constants import au as _au
from . import Sun as _Sun


gm = _Constant(
    abbrev='gm_pluto',
    name='Gravitational constant times the mass of Pluto',
    value=869.6e9,
    unit='m3 / s2',
    uncertainty=1.8e9,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

gm_system = _Constant(
    abbrev='gm_pluto_system',
    name='Gravitational constant times the mass of the Pluto system',
    value=975.5e9,
    unit='m3 / s2',
    uncertainty=1.5e9,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

mass = _Constant(
    abbrev='mass_pluto',
    name='Mass of Pluto',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_pluto and G.')

mean_radius = _Constant(
    abbrev='mean_radius_pluto',
    name='Mean radius of Pluto',
    value=1188.3e3,
    unit='m',
    uncertainty=1.6e3,
    reference='Nimmo, F., Umurhan, O., Lisse, C. M., Bierson, C. J., Lauer, '
    'T. R., Buie, M. W., Throop, H. B., Kammer, J. A., Roberts, J. H., '
    'McKinnon, W. B., Zangari, A. M., Moore, J. M., Stern, S. A., Young, '
    'L. A., Weaver, H. A., Olkin, C. B., & Ennico, K. (2017). Mean radius and '
    'shape of Pluto and Charon from New Horizons images. Icarus, 287, 12–29. '
    'https://doi.org/10.1016/j.icarus.2016.06.027')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_pluto',
    name='Volume equivalent radius of Pluto',
    value=mean_radius.value,
    unit='m',
    uncertainty=mean_radius.uncertainty,
    reference='Equal to mean_radius_pluto')

volume = _Constant(
    abbrev='volume_pluto',
    name='Volume of Pluto',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m3',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_pluto')

mean_density = _Constant(
    abbrev='mean_density_pluto',
    name='Mean density of Pluto',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_pluto and volume_equivalent_radius_pluto.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_pluto',
    name='Surface gravity of Pluto, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_pluto and mean_radius_pluto.')

angular_velocity = _Constant(
    abbrev='angular_velocity_pluto',
    name='Angular spin rate of Pluto',
    value=1.1385591834674098e-05,
    unit='rad / s',
    uncertainty=0.,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

rotational_period = _Constant(
    abbrev='rotational_period_pluto',
    name='Rotational period of Pluto',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_pluto')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_pluto',
    name='Semimajor axis of the orbit of Pluto about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=39.48211675,
    unit='au',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_pluto',
    name='Eccentricity of the orbit of Pluto about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=0.24882730,
    unit='',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_pluto',
    name='Inclination of the orbit of Pluto about the Sun, with respect to '
    'the mean ecliptic and equinox of J2000',
    value=17.14001206,
    unit='degrees',
    uncertainty=0.,
    reference='Standish, M., & Williams, J. (2012). Orbital Ephemerides of '
    'the Sun, Moon, and Planets. In Explanatory Supplement to the '
    'Astronomical Almanac (Third edition, pp. 305–342). University Science '
    'Books.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_pluto',
    name='Orbital angular velocity of the Pluto system',
    value=_np.sqrt((_Sun.gm.value + gm_system.value) /
                   (_au.value * orbit_semimajor_axis.value)**3),
    unit='rad / s',
    uncertainty=_np.sqrt(
        _Sun.gm.uncertainty**2 / 4. / (_Sun.gm.value + gm_system.value) /
        (_au.value * orbit_semimajor_axis.value)**3 +
        gm_system.uncertainty**2 / 4. / (_Sun.gm.value + gm_system.value) /
        (_au.value * orbit_semimajor_axis.value)**3 +
        9. * (_au.value * orbit_semimajor_axis.uncertainty)**2 *
        (_Sun.gm.value + gm_system.value) / 4. /
        (_au.value * orbit_semimajor_axis.value)**5),
    reference="Approximated using Kepler's third law, gm_sun, gm_pluto_system "
    'and orbit_semimajor_axis_pluto')

orbit_period = _Constant(
    abbrev='orbit_period_pluto',
    name='Orbital period of the Pluto system',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_pluto')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'rotational_period', 'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_period']
