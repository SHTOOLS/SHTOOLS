"""
pyshtools constants for the asteroid (16) Psyche.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_psyche',
    name='Gravitational constant times the mass of (16) Psyche',
    value=1.601e9,
    unit='m3 / s2',
    uncertainty=0.017e9,
    reference='Farnocchia, D., Fuentes-Muñoz, O., Park, R. S., Baer, J., and '
    'Chesley, S. R. (2024). Mass, Density, and Radius of Asteroid (16) Psyche '
    'from High-precision Astrometry. The Astronomical Journal, 168(1), 21. '
    'https://doi.org/10.3847/1538-3881/ad50ca')

mass = _Constant(
    abbrev='mass_psyche',
    name='Mass of (16) Psyche',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_psyche and G.')


mean_radius = _Constant(
    abbrev='mean_radius_psyche',
    name='Mean radius of (16) Psyche',
    value=109480.4,
    unit='m',
    uncertainty=0.,
    reference='Shepard2021_shape: Wieczorek, M. (2024). Spherical harmonic '
    'models of the shape of asteroid (16) Psyche (1.1) [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.12522382')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_psyche',
    name='Volume equivalent radius of (16) Psyche',
    value=111380.0,
    unit='m',
    uncertainty=0.,
    reference='Computed using Shepard2021_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_psyche',
    name='Volume of (16) Psyche',
    value=5787758983575747.0,
    unit='m3',
    uncertainty=0.,
    reference='Computed using Shepard2021_shape and SHCoeffs.volume()')

mean_density = _Constant(
    abbrev='mean_density_psyche',
    name='Mean density of (16) Psyche',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_psyche and volume_equivalent_radius_psyche.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_psyche',
    name='Gravity at the mean radius of (16) Psyche, ignoring rotation',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_psyche and mean_radius_psyche.')

angular_velocity = _Constant(
    abbrev='angular_velocity_psyche',
    name='Angular spin rate of (16) Psyche',
    value=2. * _np.pi / (4.195948 * 60. * 60.),
    unit='rad / s',
    uncertainty=2. * _np.pi / (60. * 60.) * 0.000001 / 4.195948**2,
    reference='Shepard, M. K., de Kleer, K., Cambioni, S., Taylor, P. A., '
    'Virkki, A. K., Rívera-Valentin, E. G., Rodriguez Sanchez-Vahamonde, C., '
    'Fernanda Zambrano-Marin, L., Magri, C., Dunham, D., Moore, J., & '
    'Camarca, M. (2021). Asteroid 16 Psyche: Shape, Features, and Global Map. '
    'The Planetary Science Journal, 2(4), 125. '
    'https://doi.org/10.3847/PSJ/abfdba')

rotational_period = _Constant(
    abbrev='rotational_period_psyche',
    name='Rotational period of (16) Psyche',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_psyche')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_psyche',
    name='Semimajor axis of the orbit of (16) Psyche about the Sun',
    value=2.923662128456563,
    unit='au',
    uncertainty=6.6729E-9,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2022-Feb-03 10:45:03')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_psyche',
    name='Eccentricity of the orbit of (16) Psyche about the Sun',
    value=0.1341931447944021,
    unit='',
    uncertainty=2.8878E-8,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2022-Feb-03 10:45:03')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_psyche',
    name='Inclination of the orbit of (16) Psyche about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=3.096815244138823,
    unit='degrees',
    uncertainty=3.6707E-6,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2022-Feb-03 10:45:03')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_psyche',
    name='Orbital angular velocity of (16) Psyche about the Sun',
    value=2 * _np.pi / (1825.95134138117 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=2 * _np.pi / (24 * 60 * 60) * 6.2513E-6 / 1825.95134138117**2,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2022-Feb-03 10:45:03')

orbit_period = _Constant(
    abbrev='orbit_period_psyche',
    name='Orbital period of (16) Psyche',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_psyche')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density',
           'angular_velocity', 'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'rotational_period',
           'orbit_period']
