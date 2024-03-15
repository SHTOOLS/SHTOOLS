"""
pyshtools constants for the asteroid (433) Eros.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_eros',
    name='Gravitational constant times the mass of (433) Eros',
    value=4.46275472004e5,
    unit='m3 / s2',
    uncertainty=2.19591611330,
    reference='JGE15A01: Miller, J. K., Konopliv, A. S., Antreasian, P. G., '
    'Bordi, J. J., Chesley, S., Helfrich, C. E., Owen, W. M., Wang, T. C., '
    'Williams, B. G., Yeomans, D. K., & Scheeres, D. J. (2002). Determination '
    'of Shape, Gravity, and Rotational State of Asteroid 433 Eros. Icarus, '
    '155(1), 3–17. https://doi.org/10.1006/icar.2001.6753')

mass = _Constant(
    abbrev='mass_eros',
    name='Mass of (433) Eros',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_eros and G.')

mean_radius = _Constant(
    abbrev='r_eros',
    name='Mean radius of (433) Eros',
    value=7315.9,
    unit='m',
    uncertainty=0.,
    reference='NLR_shape: Wieczorek, M. (2024). Spherical harmonic models of '
    'the shape of (433) Eros (1.0.0) [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.10820812')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='r_volume_eros',
    name='Volume equivalent radius of (433) Eros',
    value=8427.8,
    unit='m',
    uncertainty=0.,
    reference='Computed using NLR_shape and SHCoeffs.volume()')

density = _Constant(
    abbrev='density_eros',
    name='Mean density of (433) Eros',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_eros and r_volume_eros.')

g0 = _Constant(
    abbrev='g0_eros',
    name='Gravity of (433) Eros at mean planetary radius, ignoring rotation',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_eros and r_eros.')

omega = _Constant(
    abbrev='omega_vesta',
    name='Angular spin rate of Vesta',
    value=1639.389232 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='JGE15A01: Miller, J. K., Konopliv, A. S., Antreasian, P. G., '
    'Bordi, J. J., Chesley, S., Helfrich, C. E., Owen, W. M., Wang, T. C., '
    'Williams, B. G., Yeomans, D. K., & Scheeres, D. J. (2002). Determination '
    'of Shape, Gravity, and Rotational State of Asteroid 433 Eros. Icarus, '
    '155(1), 3–17. https://doi.org/10.1006/icar.2001.6753')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'gm', 'density', 'omega']
