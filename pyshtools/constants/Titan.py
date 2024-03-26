"""
pyshtools constants for Saturn's moon Titan.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_titan',
    name='Gravitational constant times the mass of Titan',
    value=8978.1383e9,
    unit='m3 / s2',
    uncertainty=0.0003e9,
    reference='Durante, D., Hemingway, D. J., Racioppa, P., Iess, L., & '
    'Stevenson, D. J. (2019). Titan’s gravity field and interior structure '
    'after Cassini. Icarus, 326, 123–132. '
    'https://doi.org/10.1016/j.icarus.2019.03.003')

mass = _Constant(
    abbrev='mass_titan',
    name='Mass of Titan',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_titan and G.')

mean_radius = _Constant(
    abbrev='mean_radius_titan',
    name='Mean radius of Titan',
    value=2574761.2,
    unit='m',
    uncertainty=17.7,
    reference='Corlies2017_shape: Corlies, P., Hayes, A. G., Birch, S. P. D., '
    'Lorenz, R., Stiles, B. W., Kirk, R., Poggiali, V., Zebker, H., & Iess, '
    'L. (2017). Titan’s Topography and Shape at the End of the Cassini '
    'Mission. Geophysical Research Letters, 44(23), 11,754-11,761. '
    'https://doi.org/10.1002/2017GL075518')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_titan',
    name='Volume equivalent radius of Titan',
    value=2574761.2,
    unit='m',
    uncertainty=0.,
    reference='Computed using Corlies2017 and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_titan',
    name='Volume of Titan',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_titan')

mean_density = _Constant(
    abbrev='mean_density_titan',
    name='Mean density of Titan',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_titan and volume_equivalent_radius_titan.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_titan',
    name='Gravity at the mean radius of Titan, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_titan and mean_radius_titan.')

angular_velocity = _Constant(
    abbrev='angular_velocity_titan',
    name='Angular spin rate of Titan',
    value=2 * _np.pi / (15.945448 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

rotational_period = _Constant(
    abbrev='rotational_period_titan',
    name='Rotational period of Titan',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_titan')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_titan',
    name='Semimajor axis of the orbit of Titan about Saturn',
    value=1221900.e3,
    unit='m',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_titan',
    name='Eccentricity of the orbit of Titan about Saturn',
    value=0.029,
    unit='',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_titan',
    name='Inclination of the orbit of Titan about Saturn with respect '
    'to the Laplace plane',
    value=0.3,
    unit='degrees',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn’s Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_titan',
    name='Orbital angular velocity of Titan about Saturn',
    value=2 * _np.pi / (15.945448 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_tilt = _Constant(
    abbrev='orbit_tilt_titan',
    name='Angle between the Titan Laplace plane and the equatorial plane '
    'of Saturn',
    value=0.6,
    unit='degrees',
    uncertainty=0.,
    reference='Jacobson, R. (2022). The Orbits of the Main Saturnian '
    'Satellites, the Saturnian System Gravity Field, and the Orientation of '
    "Saturn's Pole. The Astronomical Journal, 164, 199. "
    'https://doi.org/10.3847/1538-3881/ac90c9')

orbit_period = _Constant(
    abbrev='orbit_period_titan',
    name='Orbital period of Titan',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_titan')


__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_tilt',
           'rotational_period', 'orbit_period']
