"""
pyshtools constants for the planet Earth.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm_egm2008 = _Constant(
    abbrev='gm_egm2008',
    name='Gravitational constant times the mass of Earth for the model '
         'EGM2008, including the atmosphere',
    value=3986004.415e+8,
    unit='m3 / s2',
    uncertainty=0.0,
    reference='Pavlis N. K., S. A. Holmes, S. C. Kenyon, and J. K. Factor '
    '(2012). The development and evaluation of the Earth Gravitational '
    'Model 2008 (EGM2008). J. Geophys. Res., 117, B04406, '
    'doi:10.1029/2011JB008916.')

mass_egm2008 = _Constant(
    abbrev='mass_egm2008',
    name='Mass of Earth for the model EGM2008, including the atmosphere',
    value=gm_egm2008.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm_egm2008.uncertainty / _G.value)**2 +
                         (gm_egm2008.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_egm2008 and G.')

omega_egm2008 = _Constant(
    abbrev='omega_egm2008',
    name='Angular spin rate of Earth for the model EGM2008',
    value=7292115.0e-11,
    unit='rad / s',
    uncertainty=0.0,
    reference='Pavlis N. K., S. A. Holmes, S. C. Kenyon, and J. K. Factor '
    '(2012). The development and evaluation of the Earth Gravitational '
    'Model 2008 (EGM2008). J. Geophys. Res., 117, B04406, '
    'doi:10.1029/2011JB008916.')

a_wgs84 = _Constant(
    abbrev='a_wgs84',
    name='Semimajor axis of the WGS84 ellipsoid',
    value=6378137.0,
    unit='m',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

f_wgs84 = _Constant(
    abbrev='f_wgs84',
    name='Flattening of the WGS84 ellipsoid',
    value=1/298.257223563,
    unit='',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

gm_wgs84 = _Constant(
    abbrev='gm_wgs84',
    name='Gravitational constant times the mass of Earth for the '
         'WGS84 geodetic reference system, including the atmosphere',
    value=3986004.418e8,
    unit='m3 / s2',
    uncertainty=0.008e8,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

mass_wgs84 = _Constant(
    abbrev='mass_wgs84',
    name='Mass of Earth for the WGS84 geodetic reference system, '
    'including the atmosphere',
    value=gm_wgs84.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm_wgs84.uncertainty / _G.value)**2 +
                         (gm_wgs84.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_wgs84 and G.')

omega_wgs84 = _Constant(
    abbrev='omega_wgs84',
    name='Angular rotation rate of Earth for the WGS84 geodetic reference '
    'system',
    value=7292115.0e-11,
    unit='rad / s',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

gma_wgs84 = _Constant(
    abbrev='gma_wgs84',
    name="Gravitational constant times the mass of Earth's atmosphere "
    'for the WGS84 geodetic reference system',
    value=3.5e8,
    unit='m3 / s2',
    uncertainty=0.1e8,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

b_wgs84 = _Constant(
    abbrev='b_wgs84',
    name='Semiminor axis of the WGS84 ellipsoid',
    value=6356752.3142,
    unit='m',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

r3_wgs84 = _Constant(
    abbrev='r3_wgs84',
    name="Radius of a sphere with volume equal to that of the WGS84 "
    'ellipsoid',
    value=6371000.7900,
    unit='m',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')

u0_wgs84 = _Constant(
    abbrev='u0_wgs84',
    name="Theoretical normal gravity potential of the WGS84 ellipsoid",
    value=62636851.7146,
    unit='m2 / s2',
    uncertainty=0.0,
    reference='National Imagery and Mapping Agency (2000). Department of '
    'Defense World Geodetic System 1984: Its Definition and Relationship '
    'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
    'Mapping Agency.')
