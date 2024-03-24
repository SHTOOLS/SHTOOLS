"""
pyshtools constants for the asteroid (4) Vesta.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_vesta',
    name='Gravitational constant times the mass of (4) Vesta',
    value=317288244969.3,
    unit='m3 / s2',
    uncertainty=4064.9,
    reference='VESTA20H: Konopliv, A.S., Asmar, S.W., Park, R.S., Bills, '
    'B.G., Centinello, F., Chamberlin, A.B., Ermakov, A., Gaskell, R.W., '
    'Rambaux, N., Raymond, C.A., Russell, C.T., Smith, D.E., Tricarico, P., '
    'Zuber, M.T. (2014). The Vesta gravity field, spin pole and rotation '
    'period, landmark positions, and ephemeris from the Dawn tracking and '
    'optical data. Icarus, 240, 103-117, doi:10.1016/j.icarus.2013.09.005.')

mass = _Constant(
    abbrev='mass_vesta',
    name='Mass of Vesta',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_vesta and G.')

mean_radius = _Constant(
    abbrev='mean_radius_vesta',
    name='Mean radius of (4) Vesta',
    value=260362.3,
    unit='m',
    uncertainty=0.0,
    reference='DLR_SPG_shape: Wieczorek, M. (2024). Spherical harmonic models '
    'of the shape of the asteroid (4) Vesta [DLR SPG] (1.0.0) [Data set]. '
    'Zenodo. https://doi.org/10.5281/zenodo.10800929')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_vesta',
    name='Volume equivalent radius of (4) Vesta',
    value=261590.3,
    unit='m',
    uncertainty=0.,
    reference='Computed using DLR_SPG_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_vesta',
    name='Volume of (4) Vesta',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_vesta')

mean_density = _Constant(
    abbrev='mean_density_vesta',
    name='Mean density of (4) Vesta',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_vesta and volume_equivalent_radius_vesta.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_vesta',
    name='Gravity at the mean radius of (4) Vesta, ignoring rotation',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_vesta and mean_radius_vesta.')

angular_velocity = _Constant(
    abbrev='angular_velocity_vesta',
    name='Angular spin rate of (4) Vesta',
    value=1617.3331279 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='VESTA20H: Konopliv, A.S., Asmar, S.W., Park, R.S., Bills, '
    'B.G., Centinello, F., Chamberlin, A.B., Ermakov, A., Gaskell, R.W., '
    'Rambaux, N., Raymond, C.A., Russell, C.T., Smith, D.E., Tricarico, P., '
    'Zuber, M.T. (2014). The Vesta gravity field, spin pole and rotation '
    'period, landmark positions, and ephemeris from the Dawn tracking and '
    'optical data. Icarus, 240, 103-117, doi:10.1016/j.icarus.2013.09.005.')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_vesta',
    name='Semimajor axis of the orbit of (4) Vesta about the Sun',
    value=2.361922083328795,
    unit='au',
    uncertainty=1.5688E-9,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:15:57')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_vesta',
    name='Eccentricity of the orbit of (4) Vesta about the Sun',
    value=0.08944909117827099,
    unit='',
    uncertainty=2.5002E-10,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:15:57')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_vesta',
    name='Inclination of the orbit of (4) Vesta about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=7.142176968213055,
    unit='degrees',
    uncertainty=2.1708E-7,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:15:57')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_vesta',
    name='Orbital angular velocity of (4) Vesta about the Sun',
    value=2 * _np.pi / (1325.857278061479 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=2 * _np.pi / (24 * 60 * 60) * 1.3209E-6 / 1325.857278061479**2,
    reference='Park, R., Folkner, W., Williams, J., & Boggs, D. (2021). The '
    'JPL Planetary and Lunar Ephemerides DE440 and DE441. The Astronomical '
    'Journal, 161, 105, https://doi.org/10.3847/1538-3881/abd414. '
    'Accessed via JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:15:57')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity']
