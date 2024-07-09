"""
pyshtools constants for Jupiter's moon Io.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_io',
    name='Gravitational constant times the mass of Io',
    value=5959.91e9,
    unit='m3 / s2',
    uncertainty=0.28e09,
    reference='Anderson, J. D., Jacobson, R. A., Lau, E. L., Moore, W. B., & '
    "Schubert, G. (2001). Io's gravity field and interior structure. J. "
    'Geophys. Res., 106, 32963–32969. https://doi.org/10.1029/2000JE001367')

mass = _Constant(
    abbrev='mass_io',
    name='Mass of Io',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_io and G.')

mean_radius = _Constant(
    abbrev='mean_radius_io',
    name='Mean radius of Io',
    value=1821.49e3,
    unit='m',
    uncertainty=500.,
    reference='Thomas, P. C., Davies, M. E., Colvin, T. R., Oberst, J., '
    'Schuster, P., Neukum, G., Carr, M. H., McEwen, A., Schubert, G., & '
    'Belton, M. J. S. (1998). The Shape of Io from Galileo Limb Measurements. '
    'Icarus, 135(1), 175–180. https://doi.org/10.1006/icar.1998.5987')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_io',
    name='Volume equivalent radius of Io',
    value=mean_radius.value,
    unit='m',
    uncertainty=mean_radius.uncertainty,
    reference='Equal to mean_radius_io')

volume = _Constant(
    abbrev='volume_io',
    name='Volume of Io',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m3',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_io')

mean_density = _Constant(
    abbrev='mean_density_io',
    name='Mean density of Io',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_io and volume_equivalent_radius_io.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_io',
    name='Surface gravity of Io, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_io and mean_radius_io.')

angular_velocity = _Constant(
    abbrev='angular_velocity_Io',
    name='Angular spin rate of Io',
    value=2 * _np.pi / (1.762732 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

rotational_period = _Constant(
    abbrev='rotational_period_io',
    name='Rotational period of Io',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_io')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_io',
    name='Semimajor axis of the orbit of Io about Jupiter',
    value=421800.e3,
    unit='m',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_io',
    name='Eccentricity of the orbit of Io about Jupiter',
    value=0.004,
    unit='',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_io',
    name='Inclination of the orbit of Io about Jupiter with respect '
    'to the Laplace plane',
    value=0.0,
    unit='degrees',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_io',
    name='Orbital angular velocity of Io about Jupiter',
    value=2 * _np.pi / (1.762732 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_tilt = _Constant(
    abbrev='orbit_tilt_io',
    name='Angle between the Io Laplace plane and the equatorial plane '
    'of Jupiter',
    value=0.0,
    unit='degrees',
    uncertainty=0.,
    reference='R. A. Jacobson (2021), The Orbits of the Regular Jovian '
    'Satellites and the Orientation of the Pole of Jupiter, personal '
    'communication to Horizons/NAIF. Accessed via JPL Solar System Dynamics, '
    'https://ssd.jpl.nasa.gov, JUP365.')

orbit_period = _Constant(
    abbrev='orbit_period_io',
    name='Orbital period of Io',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    orbit_angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_io')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_tilt',
           'rotational_period', 'orbit_period']
