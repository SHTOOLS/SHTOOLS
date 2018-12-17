"""
pyshtools constants for the planet Mars.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm_mars = _Constant(
    abbrev='gm_mars',
    name='Gravitational constant times the mass of Mars',
    value=0.4282837581575610e+14,
    unit='m3 / s2',
    uncertainty=0.18167460e+6,
    reference='Konopliv A. S., R. S. Park, W. M. Folkner (2016). '
    'An improved JPL Mars gravity field and orientation from Mars orbiter '
    'and lander tracking data, Icarus, 274, 253-260, '
    'doi:10.1016/j.icarus.2016.02.052')

mass_mars = _Constant(
    abbrev='mass_mars',
    name='Mass of Mars',
    value=gm_mars.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm_mars.uncertainty / _G.value)**2 +
                         (gm_mars.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_mars and G.')

r_mars = _Constant(
    abbrev='r_mars',
    name='Mean radius of Mars',
    value=3389.500e3,
    unit='m',
    uncertainty=0.0,
    reference='MarsTopo2600: Wieczorek, M. A. (2015). Gravity and '
    'topography of the terrestrial planets. In T. Spohn & G. Schubert '
    '(Eds.), Treatise on Geophysics, 2nd ed., Vol. 10, pp. 153-193). '
    'Oxford, Elsevier-Pergamon, doi:10.1016/B978-0-444-53802-4.00169-X.')

density_mars = _Constant(
    abbrev='density_mars',
    name='Mean density of Mars',
    value=3 * mass_mars.value / (_np.pi * 4 * r_mars.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass_mars.uncertainty /
                         (_np.pi * 4 * r_mars.value**3))**2
                         + (3 * 3 * mass_mars.value *
                         r_mars.uncertainty /
                         (_np.pi * 4 * r_mars.value**4))**2
                         ),
    reference='Derived from mass_mars and r_mars.')

g0_mars = _Constant(
    abbrev='g0_mars',
    name='Mean surface gravity of Mars at mean planetary radius, '
    'ignoring rotation and tides',
    value=gm_mars.value / r_mars.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm_mars.uncertainty / r_mars.value**2)**2
                         + (2 * gm_mars.value * r_mars.uncertainty
                         / r_mars.value**3)**2
                         ),
    reference='Derived from gm_mars and r_mars.')

omega_mars = _Constant(
    abbrev='omega_mars',
    name='Angular spin rate of Mars',
    value=350.891985307 * 2 * _np.pi / 360 / (24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.000000003 * 2 * _np.pi / 360 / (24 * 60 * 60),
    reference='Konopliv A. S., R. S. Park, W. M. Folkner (2016). '
    'An improved JPL Mars gravity field and orientation from Mars orbiter '
    'and lander tracking data, Icarus, 274, 253-260, '
    'doi:10.1016/j.icarus.2016.02.052')

a_mars = _Constant(
    abbrev='a_mars',
    name='Semimajor axis of the Mars reference ellipsoid',
    value=3395428.0,
    unit='m',
    uncertainty=19.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

b_mars = _Constant(
    abbrev='b_mars',
    name='Semiminor axis of the Mars reference ellipsoid',
    value=3377678.0,
    unit='m',
    uncertainty=19.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

f_mars = _Constant(
    abbrev='f_mars',
    name='Flattening of the Mars reference ellipsoid',
    value=(a_mars.value - b_mars.value) / a_mars.value,
    unit='',
    uncertainty=_np.sqrt((a_mars.uncertainty * (a_mars.value - b_mars.value)
                         / a_mars.value**2)**2
                         + (_np.sqrt(a_mars.uncertainty**2 +
                            b_mars.uncertainty**2) / a_mars.value)**2
                         ),
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')

u0_mars = _Constant(
    abbrev='u0_mars',
    name='Theoretical normal gravity potential of the reference ellipsoid',
    value=12654875.0,
    unit='m2 / s2',
    uncertainty=69.0,
    reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
    'A new reference equipotential surface, and reference ellipsoid for '
    'the planet Mars. Earth, Moon, and Planets, 106, 1-13, '
    'doi:10.1007/s11038-009-9342-7.')
