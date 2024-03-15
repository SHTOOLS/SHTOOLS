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
    abbrev='r_callisto',
    name='Mean radius of Callisto',
    value=2410.3e3,
    unit='m',
    uncertainty=1.5e3,
    reference='Anderson, J. D., Jacobson, R. A., McElrath, T. P., Moore, '
    'W. B., & Schubert, G. (2001). Shape, mean radius, gravity field, and '
    'interior structure of Callisto. Icarus, 153(1), 157–161. '
    'https://doi.org/10.1006/icar.2001.6664')

r = mean_radius

density = _Constant(
    abbrev='density_callisto',
    name='Mean density of Callisto',
    value=3 * mass.value / (_np.pi * 4 * mean_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * mean_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         mean_radius.uncertainty /
                         (_np.pi * 4 * mean_radius.value**4))**2
                         ),
    reference='Derived from mass_callisto and r_callisto.')

g0 = _Constant(
    abbrev='g0_callisto',
    name='Surface gravity of Callisto, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_callisto and r_callisto.')

omega = _Constant(
    abbrev='omega_callisto',
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


__all__ = ['gm', 'mass', 'mean_radius', 'r', 'gm', 'density', 'omega']
