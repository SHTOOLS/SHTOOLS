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
    abbrev='r_ganymede',
    name='Mean radius of Ganymede',
    value=2632630.,
    unit='m',
    uncertainty=100,
    reference='Zubarev, A., Nadezhdina, I., Oberst, J., Hussmann, H., & '
    'Stark, A. (2015). New Ganymede control point network and global shape '
    'model. Planetary and Space Science, 117, 246–249. '
    'https://doi.org/10.1016/j.pss.2015.06.022')

r = mean_radius

density = _Constant(
    abbrev='density_ganymede',
    name='Mean density of Ganymede',
    value=3 * mass.value / (_np.pi * 4 * mean_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * mean_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         mean_radius.uncertainty /
                         (_np.pi * 4 * mean_radius.value**4))**2
                         ),
    reference='Derived from mass_ganymede and r_ganymede.')

g0 = _Constant(
    abbrev='g0_ganymede',
    name='Surface gravity of Ganymede, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_ganymede and r_ganymede.')

omega = _Constant(
    abbrev='omega_Ganymede',
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


__all__ = ['gm', 'mass', 'mean_radius', 'r', 'gm', 'density', 'omega']
