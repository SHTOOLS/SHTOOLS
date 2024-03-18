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
    unit='m',
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

omega = _Constant(
    abbrev='omega_enceladus',
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

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'omega']
