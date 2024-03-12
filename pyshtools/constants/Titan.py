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
    abbrev='r_titan',
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

density = _Constant(
    abbrev='density_titan',
    name='Mean density of Titan',
    value=3 * mass.value / (_np.pi * 4 * mean_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * mean_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         mean_radius.uncertainty /
                         (_np.pi * 4 * mean_radius.value**4))**2
                         ),
    reference='Derived from mass_titan and r_titan.')

g0 = _Constant(
    abbrev='g0_titan',
    name='Surface gravity of Titan, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_titan and r_titan.')

omega = _Constant(
    abbrev='omega_titan',
    name='Angular spin rate of Titan',
    value=22.5769768 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='Archinal, B. A., Acton, C. H., A’Hearn, M. F., Conrad, A., '
    'Consolmagno, G. J., Duxbury, T., Hestroffer, D., Hilton, J. L., Kirk, '
    'R. L., Klioner, S. A., McCarthy, D., Meech, K., Oberst, J., Ping, J., '
    'Seidelmann, P. K., Tholen, D. J., Thomas, P. C., & Williams, I. P. '
    '(2018). Report of the IAU Working Group on Cartographic Coordinates and '
    'Rotational Elements: 2015. Celestial Mechanics and Dynamical Astronomy, '
    '130(3), 22. https://doi.org/10.1007/s10569-017-9805-5')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'gm', 'density', 'omega']
