"""
pyshtools constants for Jupiter's moon Ganymede.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_ganymede',
    name='Gravitational constant times the mass of Ganymede',
    value=9.8878041807018262e+12,
    unit='m3 / s2',
    uncertainty=0.10256634173857643e+09,
    reference='Gomez Casajus, L., Ermakov, A. I., Zannoni, M., Keane, J. T., '
    'Stevenson, D., Buccino, D. R., Durante, D., Parisi, M., Park, R. S., '
    'Tortora, P., Bolton, S. J. (2022). Gravity Field of Ganymede After the '
    'Juno Extended Mission. Geophysical Research Letters, 49(24), '
    'e2022GL099475, doi:10.1029/2022GL099475.')

mass = _Constant(
    abbrev='mass_ganymede',
    name='Mass of Ganymede',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_ganymede and G.')

mean_radius = _Constant(
    abbrev='mean_radius_ganymede',
    name='Mean radius of Ganymede',
    value=2632630.,
    unit='m',
    uncertainty=100,
    reference='Zubarev, A., Nadezhdina, I., Oberst, J., Hussmann, H., & '
    'Stark, A. (2015). New Ganymede control point network and global shape '
    'model. Planetary and Space Science, 117, 246–249. '
    'https://doi.org/10.1016/j.pss.2015.06.022')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_ganymede',
    name='Volume equivalent radius of Ganymede',
    value=mean_radius.value,
    unit='m',
    uncertainty=mean_radius.uncertainty,
    reference='Equal to mean_radius_ganymede')

volume = _Constant(
    abbrev='volume_ganymede',
    name='Volume of Ganymede',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_ganymede')

mean_density = _Constant(
    abbrev='mean_density_ganymede',
    name='Mean density of Ganymede',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_ganymede and '
    'volume_equivalent_radius_ganymede.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_ganymede',
    name='Surface gravity of Ganymede, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_ganymede and mean_radius_ganymede.')

angular_velocity = _Constant(
    abbrev='angular_velocity_Ganymede',
    name='Angular spin rate of Ganymede',
    value=50.3176081 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='Archinal, B. A., Acton, C. H., A’Hearn, M. F., Conrad, A., '
    'Consolmagno, G. J., Duxbury, T., Hestroffer, D., Hilton, J. L., Kirk, '
    'R. L., Klioner, S. A., McCarthy, D., Meech, K., Oberst, J., Ping, J., '
    'Seidelmann, P. K., Tholen, D. J., Thomas, P. C., & Williams, I. P. '
    '(2018). Report of the IAU Working Group on Cartographic Coordinates and '
    'Rotational Elements: 2015. Celestial Mechanics and Dynamical Astronomy, '
    '130(3), 22. https://doi.org/10.1007/s10569-017-9805-5')

rotational_period = _Constant(
    abbrev='rotational_period_ganymede',
    name='Rotational period of Ganymede',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_ganymede')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_ganymede',
    name='Semimajor axis of the orbit of Ganymede about Jupiter',
    value=1070400.e3,
    unit='m',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_ganymede',
    name='Eccentricity of the orbit of Ganymede about Jupiter',
    value=0.001,
    unit='',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_ganymede',
    name='Inclination of the orbit of Ganymede about Jupiter with respect '
    'to the Laplace plane',
    value=0.2,
    unit='degrees',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_ganymede',
    name='Orbital angular velocity of Ganymede about Jupiter',
    value=2 * _np.pi / (7.155588 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_tilt = _Constant(
    abbrev='orbit_tilt_ganymede',
    name='Angle between the Ganymede Laplace plane and the equatorial plane '
    'of Jupiter',
    value=0.1,
    unit='degrees',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_period = _Constant(
    abbrev='orbit_period_ganymede',
    name='Orbital period of Ganymede',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_ganymede')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_tilt',
           'rotational_period', 'orbit_period']
