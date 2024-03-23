"""
pyshtools constants for the asteroid (1) Ceres.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_ceres',
    name='Gravitational constant times the mass of (1) Ceres',
    value=62629053612.1,
    unit='m3 / s2',
    uncertainty=350000.0,
    reference='CERES18D: Konopliv, A.S., Park, R.S., Vaughan, A.T., Bills, '
    'B.G., Asmar, S.W., Ermakov, A.I., Rambaux, N., Raymond, C.A., '
    'Castillo-Rogez, J.C., Russell, C.T., Smith, D.E., Zuber, M.T. (2018). '
    'The Ceres gravity field, spin pole, rotation period and orbit from the '
    'Dawn radiometric tracking and optical data, Icarus, 299, 411-429, '
    'doi:10.1016/j.icarus.2017.08.005.')

mass = _Constant(
    abbrev='mass_ceres',
    name='Mass of Ceres',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_ceres and G.')

mean_radius = _Constant(
    abbrev='mean_radius_ceres',
    name='Mean radius of (1) Ceres',
    value=469461.8,
    unit='m',
    uncertainty=0.,
    reference='JPL_SPC_shape: Wieczorek, M. (2024). Spherical harmonic models '
    'of the shape of asteroid (1) Ceres [JPL SPC] (1.0.0) [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.10812848')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_ceres',
    name='Volume equivalent radius of (1) Ceres',
    value=469725.0,
    unit='m',
    uncertainty=0.,
    reference='Computed using JPL_SPC_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_ceres',
    name='Volume of (1) Ceres',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_ceres')

mean_density = _Constant(
    abbrev='mean_density_ceres',
    name='Mean density of (1) Ceres',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_ceres and volume_equivalent_radius_ceres.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_ceres',
    name='Gravity at the mean radius of (1) Ceres, ignoring rotation',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_ceres and mean_radius_ceres.')

angular_velocity = _Constant(
    abbrev='angular_velocity_ceres',
    name='Angular spin rate of (1) Ceres',
    value=952.1532635 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.000002 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    reference='CERES18D: Konopliv, A.S., Park, R.S., Vaughan, A.T., Bills, '
    'B.G., Asmar, S.W., Ermakov, A.I., Rambaux, N., Raymond, C.A., '
    'Castillo-Rogez, J.C., Russell, C.T., Smith, D.E., Zuber, M.T. (2018). '
    'The Ceres gravity field, spin pole, rotation period and orbit from the '
    'Dawn radiometric tracking and optical data, Icarus, 299, 411-429, '
    'doi:10.1016/j.icarus.2017.08.005.')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_ceres',
    name='Semimajor axis of the orbit of (1) Ceres about the Sun',
    value=2.767254360873952,
    unit='au',
    uncertainty=1.0197E-11,
    reference='JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:04:44')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_ceres',
    name='Eccentricity of the orbit of (1) Ceres about the Sun',
    value=0.0789125317658808,
    unit='',
    uncertainty=4.7845E-12,
    reference='JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:04:44')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_ceres',
    name='Inclination of the orbit of (1) Ceres about the Sun, with respect '
    'to the mean ecliptic and equinox of J2000',
    value=10.5868796009696,
    unit='degrees',
    uncertainty=4.6112E-9,
    reference='JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:04:44')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_ceres',
    name='Orbital angular velocity of (1) Ceres about the Sun',
    value=2 * _np.pi / (1681.403889244621 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=2 * _np.pi / (24 * 60 * 60) * 9.2941E-9 / 1681.403889244621**2,
    reference='JPL Solar System Dynamics, https://ssd.jpl.nasa.gov, '
    'solution date: 2021-Apr-13 11:04:44')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density',
           'angular_velocity', 'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity']
