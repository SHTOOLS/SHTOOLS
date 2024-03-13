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
    name='Gravitational constant times the mass of Ceres',
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
    abbrev='r_ceres',
    name='Mean radius of Ceres',
    value=469461.8,
    unit='m',
    uncertainty=0.,
    reference='JPL_SPC_shape: Wieczorek, M. (2024). Spherical harmonic models '
    'of the shape of asteroid (1) Ceres [JPL SPC] (1.0.0) [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.10812848')

r = mean_radius

density = _Constant(
    abbrev='density_ceres',
    name='Mean density of Ceres',
    value=3 * mass.value / (_np.pi * 4 * mean_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * mean_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         mean_radius.uncertainty /
                         (_np.pi * 4 * mean_radius.value**4))**2
                         ),
    reference='Derived from mass_ceres and r_ceres.')

g0 = _Constant(
    abbrev='g0_ceres',
    name='Surface gravity of Ceres, ignoring rotation',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_ceres and r_ceres.')

omega = _Constant(
    abbrev='omega_ceres',
    name='Angular spin rate of Ceres',
    value=952.1532635 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.000002 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    reference='CERES18D: Konopliv, A.S., Park, R.S., Vaughan, A.T., Bills, '
    'B.G., Asmar, S.W., Ermakov, A.I., Rambaux, N., Raymond, C.A., '
    'Castillo-Rogez, J.C., Russell, C.T., Smith, D.E., Zuber, M.T. (2018). '
    'The Ceres gravity field, spin pole, rotation period and orbit from the '
    'Dawn radiometric tracking and optical data, Icarus, 299, 411-429, '
    'doi:10.1016/j.icarus.2017.08.005.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'gm', 'density', 'omega']
