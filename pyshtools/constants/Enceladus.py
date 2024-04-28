"""
pyshtools constants for Saturn's moon Enceladus.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_enceladus',
    name='Gravitational constant times the mass of Enceladus',
    value=7.210443e9,
    unit='m3 / s2',
    uncertainty=0.000030e9,
    reference='Park, R. S., Mastrodemos, N., Jacobson, R. A., Berne, A., '
    'Vaughan, A. T., Hemingway, D. J., Leonard, E. J., Castillo-Rogez, J. C., '
    'Cockell, C. S., Keane, J. T., Konopliv, A. S., Nimmo, F., Riedel, J. E., '
    'Simons, M., & Vance, S. (2024). The Global Shape, Gravity Field, and '
    'Libration of Enceladus. Journal of Geophysical Research: Planets, '
    '129(1), e2023JE008054. https://doi.org/10.1029/2023JE008054')

mass = _Constant(
    abbrev='mass_enceladus',
    name='Mass of Enceladus',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_enceladus and G.')

mean_radius = _Constant(
    abbrev='r_enceladus',
    name='Mean radius of Enceladus',
    value=251967.3,
    unit='m',
    uncertainty=0.,
    reference='Park, R. S., Mastrodemos, N., Jacobson, R. A., Berne, A., '
    'Vaughan, A. T., Hemingway, D. J., Leonard, E. J., Castillo-Rogez, J. C., '
    'Cockell, C. S., Keane, J. T., Konopliv, A. S., Nimmo, F., Riedel, J. E., '
    'Simons, M., & Vance, S. (2024). The Global Shape, Gravity Field, and '
    'Libration of Enceladus. Journal of Geophysical Research: Planets, '
    '129(1), e2023JE008054. https://doi.org/10.1029/2023JE008054')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_enceladus',
    name='Volume equivalent radius of Enceladus',
    value=251985.3,
    unit='m',
    uncertainty=0.,
    reference='Computed using JPL_SPC_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_enceladus',
    name='Volume of Enceladus',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m3',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_enceladus')

mean_density = _Constant(
    abbrev='mean_density_enceladus',
    name='Mean density of Enceladus',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_enceladus and '
    'volume_equivalent_radius_enceladus.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_enceladus',
    name='Gravity at the mean radius of Enceladus, ignoring rotation and '
    'tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_enceladus and mean_radius_enceladus.')

angular_velocity = _Constant(
    abbrev='angular_velocity_enceladus',
    name='Angular spin rate of Enceladus',
    value=262.7318870466 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='Park, R. S., Mastrodemos, N., Jacobson, R. A., Berne, A., '
    'Vaughan, A. T., Hemingway, D. J., Leonard, E. J., Castillo-Rogez, J. C., '
    'Cockell, C. S., Keane, J. T., Konopliv, A. S., Nimmo, F., Riedel, J. E., '
    'Simons, M., & Vance, S. (2024). The Global Shape, Gravity Field, and '
    'Libration of Enceladus. Journal of Geophysical Research: Planets, '
    '129(1), e2023JE008054. https://doi.org/10.1029/2023JE008054')

rotational_period = _Constant(
    abbrev='rotational_period_enceladus',
    name='Rotational period of Enceladus',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_enceladus')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_enceladus',
    name='Semimajor axis of the orbit of Enceladus about Saturn',
    value=238400e3,
    unit='m',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_enceladus',
    name='Eccentricity of the orbit of Enceladus about Saturn',
    value=0.005,
    unit='',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_enceladus',
    name='Inclination of the orbit of Enceladus about Saturn with respect '
    'to the Laplace plane',
    value=0.0,
    unit='degrees',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_enceladus',
    name='Orbital angular velocity of enceladus about Saturn',
    value=2 * _np.pi / (1.370218 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_tilt = _Constant(
    abbrev='orbit_tilt_enceladus',
    name='Angle between the Enceladus Laplace plane and the equatorial plane '
    'of Saturn',
    value=0.0,
    unit='degrees',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_period = _Constant(
    abbrev='orbit_period_enceladus',
    name='Orbital period of Enceladus',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_enceladus')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_tilt',
           'rotational_period', 'orbit_period']
