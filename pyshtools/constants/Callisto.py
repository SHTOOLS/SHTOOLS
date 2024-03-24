"""
pyshtools constants for Jupiter's moon Callisto.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_callisto',
    name='Gravitational constant times the mass of Callisto',
    value=7179.292e9,
    unit='m3 / s2',
    uncertainty=0.009e9,
    reference='Anderson, J. D., Jacobson, R. A., McElrath, T. P., Moore, '
    'W. B., & Schubert, G. (2001). Shape, mean radius, gravity field, and '
    'interior structure of Callisto. Icarus, 153(1), 157–161. '
    'https://doi.org/10.1006/icar.2001.6664')

mass = _Constant(
    abbrev='mass_callisto',
    name='Mass of Callisto',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_callisto and G.')

mean_radius = _Constant(
    abbrev='mean_radius_callisto',
    name='Mean radius of Callisto',
    value=2410.3e3,
    unit='m',
    uncertainty=1.5e3,
    reference='Anderson, J. D., Jacobson, R. A., McElrath, T. P., Moore, '
    'W. B., & Schubert, G. (2001). Shape, mean radius, gravity field, and '
    'interior structure of Callisto. Icarus, 153(1), 157–161. '
    'https://doi.org/10.1006/icar.2001.6664')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_callisto',
    name='Volume equivalent radius of Callisto',
    value=mean_radius.value,
    unit='m',
    uncertainty=mean_radius.uncertainty,
    reference='Equal to mean_radius_callisto')

volume = _Constant(
    abbrev='volume_callisto',
    name='Volume of Callisto',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_callisto')

mean_density = _Constant(
    abbrev='mean_density_callisto',
    name='Mean density of Callisto',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_callisto and '
    'volume_equivalent_radius_callisto.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_callisto',
    name='Gravity at the mean radius of Callisto, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_callisto and mean_radius_callisto.')

angular_velocity = _Constant(
    abbrev='angular_velocity_callisto',
    name='Angular spin rate of Callisto',
    value=21.5710715 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='Archinal, B. A., Acton, C. H., A’Hearn, M. F., Conrad, A., '
    'Consolmagno, G. J., Duxbury, T., Hestroffer, D., Hilton, J. L., Kirk, '
    'R. L., Klioner, S. A., McCarthy, D., Meech, K., Oberst, J., Ping, J., '
    'Seidelmann, P. K., Tholen, D. J., Thomas, P. C., & Williams, I. P. '
    '(2018). Report of the IAU Working Group on Cartographic Coordinates and '
    'Rotational Elements: 2015. Celestial Mechanics and Dynamical Astronomy, '
    '130(3), 22. https://doi.org/10.1007/s10569-017-9805-5')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_callisto',
    name='Semimajor axis of the orbit of Callisto about Jupiter',
    value=1882700.e3,
    unit='m',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_callisto',
    name='Eccentricity of the orbit of Callisto about Jupiter',
    value=0.007,
    unit='',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_callisto',
    name='Inclination of the orbit of Callisto about Jupiter with respect '
    'to the Laplace plane',
    value=0.3,
    unit='degrees',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_callisto',
    name='Orbital angular velocity of Callisto about Jupiter',
    value=2 * _np.pi / (16.690440 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_tilt = _Constant(
    abbrev='orbit_tilt_callisto',
    name='Angle between the Callisto Laplace plane and the equatorial plane '
    'of Jupiter',
    value=0.4,
    unit='degrees',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'gravity_mean_radius', 'mean_density', 'volume', 'angular_velocity',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_tilt']
