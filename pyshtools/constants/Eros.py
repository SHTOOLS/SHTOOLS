"""
pyshtools constants for the asteroid (433) Eros.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_eros',
    name='Gravitational constant times the mass of (433) Eros',
    value=4.46275472004e5,
    unit='m3 / s2',
    uncertainty=2.19591611330,
    reference='JGE15A01: Miller, J. K., Konopliv, A. S., Antreasian, P. G., '
    'Bordi, J. J., Chesley, S., Helfrich, C. E., Owen, W. M., Wang, T. C., '
    'Williams, B. G., Yeomans, D. K., & Scheeres, D. J. (2002). Determination '
    'of Shape, Gravity, and Rotational State of Asteroid 433 Eros. Icarus, '
    '155(1), 3–17. https://doi.org/10.1006/icar.2001.6753')

mass = _Constant(
    abbrev='mass_eros',
    name='Mass of (433) Eros',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_eros and G.')

mean_radius = _Constant(
    abbrev='mean_radius_eros',
    name='Mean radius of (433) Eros',
    value=7315.9,
    unit='m',
    uncertainty=0.,
    reference='NLR_shape: Wieczorek, M. (2024). Spherical harmonic models of '
    'the shape of (433) Eros (1.0.0) [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.10820812')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_eros',
    name='Volume equivalent radius of (433) Eros',
    value=8427.8,
    unit='m',
    uncertainty=0.,
    reference='Computed using NLR_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_eros',
    name='Volume of (433) Eros',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m3',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_eros')

mean_density = _Constant(
    abbrev='mean_density_eros',
    name='Mean density of (433) Eros',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_eros and volume_equivalent_radius_eros.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_eros',
    name='Gravity at the mean radius of (433) Eros, ignoring rotation',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_eros and mean_radius_eros.')

angular_velocity = _Constant(
    abbrev='angular_velocity_eros',
    name='Angular spin rate of (433) Eros',
    value=1639.389232 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='JGE15A01: Miller, J. K., Konopliv, A. S., Antreasian, P. G., '
    'Bordi, J. J., Chesley, S., Helfrich, C. E., Owen, W. M., Wang, T. C., '
    'Williams, B. G., Yeomans, D. K., & Scheeres, D. J. (2002). Determination '
    'of Shape, Gravity, and Rotational State of Asteroid 433 Eros. Icarus, '
    '155(1), 3–17. https://doi.org/10.1006/icar.2001.6753')

rotational_period = _Constant(
    abbrev='rotational_period_eros',
    name='Rotational period of (433) Eros',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_eros')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_eros',
    name='Semimajor axis of the orbit of (433) Eros about the Sun',
    value=1.458117412303767,
    unit='au',
    uncertainty=1.5706E-10,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-May-24 17:55:05')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_eros',
    name='Eccentricity of the orbit of (433) Eros about the Sun',
    value=0.2227966940876033,
    unit='',
    uncertainty=9.3813E-9,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-May-24 17:55:05')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_eros',
    name='Inclination of the orbit of (433) Eros about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=10.82792727465937,
    unit='degrees',
    uncertainty=1.158E-6,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-May-24 17:55:05')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_eros',
    name='Orbital angular velocity of (433) Eros about the Sun',
    value=2 * _np.pi / (643.1128260948039 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=2 * _np.pi / (24 * 60 * 60) * 1.0391E-7 / 643.1128260948039**2,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-May-24 17:55:05')

orbit_period = _Constant(
    abbrev='orbit_period_eros',
    name='Orbital period of (433) Eros',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_eros')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'rotational_period',
           'orbit_period']
