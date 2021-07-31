"""
pyshtools constants for the planet Mars.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_mars',
    name='Gravitational constant times the mass of Mars',
    value=0.4282837581575610e+14,
    unit='m3 / s2',
    uncertainty=0.18167460e+6,
    reference='Konopliv A. S., R. S. Park, W. M. Folkner (2016). '
    'An improved JPL Mars gravity field and orientation from Mars orbiter '
    'and lander tracking data, Icarus, 274, 253-260, '
    'doi:10.1016/j.icarus.2016.02.052')

mass = _Constant(
    abbrev='mass_mars',
    name='Mass of Mars',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_mars and G.')

mean_radius = _Constant(
    abbrev='r_mars',
    name='Mean radius of Mars',
    value=3389.500e3,
    unit='m',
    uncertainty=0.0,
    reference='MarsTopo2600: Wieczorek, M. A. (2015). Gravity and '
    'topography of the terrestrial planets. In T. Spohn & G. Schubert '
    '(Eds.), Treatise on Geophysics, 2nd ed., Vol. 10, pp. 153-193). '
    'Oxford, Elsevier-Pergamon, doi:10.1016/B978-0-444-53802-4.00169-X.')

r = mean_radius

density = _Constant(
    abbrev='density_mars',
    name='Mean density of Mars',
    value=3 * mass.value / (_np.pi * 4 * mean_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * mean_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         mean_radius.uncertainty /
                         (_np.pi * 4 * mean_radius.value**4))**2
                         ),
    reference='Derived from mass_mars and r_mars.')

g0 = _Constant(
    abbrev='g0_mars',
    name='Mean surface gravity of Mars at mean planetary radius, '
    'ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_mars and r_mars.')

omega = _Constant(
    abbrev='omega_mars',
    name='Angular spin rate of Mars',
    value=350.891985307 * 2 * _np.pi / 360 / (24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.000000003 * 2 * _np.pi / 360 / (24 * 60 * 60),
    reference='Konopliv A. S., R. S. Park, W. M. Folkner (2016). '
    'An improved JPL Mars gravity field and orientation from Mars orbiter '
    'and lander tracking data, Icarus, 274, 253-260, '
    'doi:10.1016/j.icarus.2016.02.052')

a = _Constant(
    abbrev='a_mars',
    name='Semimajor axis of the Mars reference ellipsoid',
    value=3395428.0,
    unit='m',
    uncertainty=19.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

b = _Constant(
    abbrev='b_mars',
    name='Semiminor axis of the Mars reference ellipsoid',
    value=3377678.0,
    unit='m',
    uncertainty=19.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

f = _Constant(
    abbrev='f_mars',
    name='Flattening of the Mars reference ellipsoid',
    value=(a.value - b.value) / a.value,
    unit='',
    uncertainty=_np.sqrt((a.uncertainty * (a.value - b.value)
                         / a.value**2)**2
                         + (_np.sqrt(a.uncertainty**2 +
                            b.uncertainty**2) / a.value)**2
                         ),
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

u0 = _Constant(
    abbrev='u0_mars',
    name='Theoretical normal gravity potential of the reference ellipsoid',
    value=12654875.0,
    unit='m2 / s2',
    uncertainty=69.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'density', 'g0', 'omega', 'a',
           'b', 'f', 'u0']
