"""
pyshtools constants.

This subpackage defines several constants used in analyzing gravity,
topography, and magnetic field data of the terrestrial planets. Each object is
an astropy Constant that possesses the attributes name, value, unit,
uncertainty, and reference.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np


try:
    from astropy.constants import Constant

    G = Constant(
        abbrev='G',
        name='Gravitational constant',
        value=6.67384e-11,
        unit='m3 / (kg s2)',
        uncertainty=0.00080e-11,
        reference='Mohr, P. J., B. N. Taylor, and D. B.Newell (2012). '
        'CODATA recommended values of the fundamental physical '
        'constants, 2010. Reviews of Modern Physics, 84, 1527–1605, '
        'doi:10.1103/RevModPhys.84.1527.')

    mu0 = Constant(
        abbrev='mu0',
        name='Magnetic constant',
        value=4 * _np.pi * 1.e-7,
        unit='T m / A',
        uncertainty=0.0,
        reference='Mohr, P. J., B. N. Taylor, and D. B.Newell (2012). '
        'CODATA recommended values of the fundamental physical constants, '
        '2010. Reviews of Modern Physics, 84, 1527–1605, '
        'doi:10.1103/RevModPhys.84.1527.')

    # == Mercury ==

    gm_mercury = Constant(
        abbrev='gm_mercury',
        name='Gravitational constant times the mass of Mercury',
        value=2.2031815411154894e+13,
        unit='m3 / s2',
        uncertainty=1.9361909444154922e+5,
        reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
        'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
        '(2014), The gravity field, orientation, and ephemeris of Mercury '
        'from MESSENGER observations after three years in orbit, J. Geophys. '
        'Res. Planets, 119, 2417–2436, doi:10.1002/2014JE004675.')

    mass_mercury = Constant(
        abbrev='mass_mercury',
        name='Mass of Mercury',
        value=gm_mercury.value / G.value,
        unit='kg',
        uncertainty=_np.sqrt((gm_mercury.uncertainty / G.value)**2 +
                             (gm_mercury.value * G.uncertainty / G.value**2)**2
                             ),
        reference='Derived from gm_mercury and G.')

    r_mercury = Constant(
        abbrev='r_mercury',
        name='Mean radius of Mercury',
        value=2439.40197456433e3,
        unit='m',
        uncertainty=0.0,
        reference='gtmes_150v05: Smith, D. E., M. T. Zuber, R. J. Phillips, '
        'S. C. Solomon, S. A. Hauck II, F. G. Lemoine, E. Mazarico, G. A. '
        'Neumann, S. J. Peale, J.-L. Margot, C. L. Johnson, M. H. Torrence, '
        'M. E. Perry, D. D. Rowlands, S. Goossens, J. W. Head, A. H. Taylor '
        '(2012). Gravity field and internal structure of Mercury from '
        'MESSENGER. Science, 336, 214-217, doi:10.1126/science.1218809.')

    density_mercury = Constant(
        abbrev='density_mercury',
        name='Mean density of Mercury',
        value=3 * mass_mercury.value / (_np.pi * 4 * r_mercury.value**3),
        unit='kg / m3',
        uncertainty=_np.sqrt((3 * mass_mercury.uncertainty /
                             (_np.pi * 4 * r_mercury.value**3))**2
                             + (3 * 3 * mass_mercury.value *
                             r_mercury.uncertainty /
                             (_np.pi * 4 * r_mercury.value**4))**2
                             ),
        reference='Derived from mass_mercury and r_mercury.')

    g0_mercury = Constant(
        abbrev='g0_mercury',
        name='Surface gravity of Mercury, ignoring rotation and tides',
        value=gm_mercury.value / r_mercury.value**2,
        unit='m / s2',
        uncertainty=_np.sqrt((gm_mercury.uncertainty / r_mercury.value**2)**2
                             + (2 * gm_mercury.value * r_mercury.uncertainty
                             / r_mercury.value**3)**2
                             ),
        reference='Derived from gm_mercury and r_mercury.')

    omega_mercury = Constant(
        abbrev='omega_mercury',
        name='Angular spin rate of Mercury',
        value=6.1385108 * 2 * _np.pi / 360 / (24 * 60 * 60),
        unit='rad / s',
        uncertainty=0.0,
        reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
        'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
        '(2014), The gravity field, orientation, and ephemeris of Mercury '
        'from MESSENGER observations after three years in orbit, J. Geophys. '
        'Res. Planets, 119, 2417–2436, doi:10.1002/2014JE004675.')

    omega_orbit_mercury = Constant(
        abbrev='omega_orbit_mercury',
        name='Angular rotation rate of Mercury about the Sun',
        value=2 * _np.pi / (87.969216879 * 24 * 60 * 60),
        unit='rad / s',
        uncertainty=6 * 2 * _np.pi / (87.969216879 * 24 * 60 * 60)**2,
        reference='ggmes_100v07: Mazarico, E., A. Genova, S. Goossens, F. G. '
        'Lemoine, G. A. Neumann, M. T. Zuber, D. E. Smith, and S. C. Solomon '
        '(2014), The gravity field, orientation, and ephemeris of Mercury '
        'from MESSENGER observations after three years in orbit, J. Geophys. '
        'Res. Planets, 119, 2417–2436, doi:10.1002/2014JE004675.')

    # == Venus ==

    gm_venus = Constant(
        abbrev='gm_venus',
        name='Gravitational constant times the mass of Venus',
        value=324858592079000.,
        unit='m3 / s2',
        uncertainty=6376000.0,
        reference='MGNP180U: Konopliv A. S., W. B. Banerdt, and W. L. Sjogren '
        '(1999) Venus gravity: 180th degree and order model. Icarus 139: 3–18.'
        'doi:10.1006/icar.1999.6086.')

    mass_venus = Constant(
        abbrev='mass_venus',
        name='Mass of Venus',
        value=gm_venus.value / G.value,
        unit='kg',
        uncertainty=_np.sqrt((gm_venus.uncertainty / G.value)**2 +
                             (gm_venus.value * G.uncertainty / G.value**2)**2
                             ),
        reference='Derived from gm_venus and G.')

    r_venus = Constant(
        abbrev='r_venus',
        name='Mean radius of Venus',
        value=6051.878e3,
        unit='m',
        uncertainty=0.0,
        reference='VenusTopo719: Wieczorek, M. A. (2015). Gravity and '
        'topography of the terrestrial planets. In T. Spohn & G. Schubert '
        '(Eds.), Treatise on Geophysics, 2nd ed., Vol. 10, pp. 153–193). '
        'Oxford, Elsevier-Pergamon, doi:10.1016/B978-0-444-53802-4.00169-X.')

    density_venus = Constant(
        abbrev='density_venus',
        name='Mean density of Venus',
        value=3 * mass_venus.value / (_np.pi * 4 * r_venus.value**3),
        unit='kg / m3',
        uncertainty=_np.sqrt((3 * mass_venus.uncertainty /
                             (_np.pi * 4 * r_venus.value**3))**2
                             + (3 * 3 * mass_venus.value *
                             r_venus.uncertainty /
                             (_np.pi * 4 * r_venus.value**4))**2
                             ),
        reference='Derived from mass_venus and r_venus.')

    g0_venus = Constant(
        abbrev='g0_venus',
        name='Surface gravity of Venus, ignoring rotation and tides',
        value=gm_venus.value / r_venus.value**2,
        unit='m / s2',
        uncertainty=_np.sqrt((gm_venus.uncertainty / r_venus.value**2)**2
                             + (2 * gm_venus.value * r_venus.uncertainty
                             / r_venus.value**3)**2
                             ),
        reference='Derived from gm_venus and r_venus.')

    omega_venus = Constant(
        abbrev='omega_venus',
        name='Angular spin rate of Venus',
        value=-2 * _np.pi / (243.0200 * 24 * 60 * 60),
        unit='rad / s',
        uncertainty=0.0002 * 2 * _np.pi / (243.0200**2 * 24 * 60 * 60),
        reference='MGNP180U: Konopliv A. S., W. B. Banerdt, and W. L. Sjogren '
        '(1999) Venus gravity: 180th degree and order model. Icarus 139: 3–18.'
        'doi:10.1006/icar.1999.6086.')

    # == Earth ==

    gm_egm2008 = Constant(
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

    mass_egm2008 = Constant(
        abbrev='mass_egm2008',
        name='Mass of Earth for the model EGM2008, including the atmosphere',
        value=gm_egm2008.value / G.value,
        unit='kg',
        uncertainty=_np.sqrt((gm_egm2008.uncertainty / G.value)**2 +
                             (gm_egm2008.value * G.uncertainty / G.value**2)**2
                             ),
        reference='Derived from gm_egm2008 and G.')

    omega_egm2008 = Constant(
        abbrev='omega_egm2008',
        name='Angular spin rate of Earth for the model EGM2008',
        value=7292115.0e-11,
        unit='rad / s',
        uncertainty=0.0,
        reference='Pavlis N. K., S. A. Holmes, S. C. Kenyon, and J. K. Factor '
        '(2012). The development and evaluation of the Earth Gravitational '
        'Model 2008 (EGM2008). J. Geophys. Res., 117, B04406, '
        'doi:10.1029/2011JB008916.')

    a_wgs84 = Constant(
        abbrev='a_wgs84',
        name='Semimajor axis of the WGS84 ellipsoid',
        value=6378137.0,
        unit='m',
        uncertainty=0.0,
        reference='National Imagery and Mapping Agency (2000). Department of '
        'Defense World Geodetic System 1984: Its Definition and Relationship '
        'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
        'Mapping Agency.')

    f_wgs84 = Constant(
        abbrev='f_wgs84',
        name='Flattening of the WGS84 ellipsoid',
        value=1/298.257223563,
        unit='',
        uncertainty=0.0,
        reference='National Imagery and Mapping Agency (2000). Department of '
        'Defense World Geodetic System 1984: Its Definition and Relationship '
        'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
        'Mapping Agency.')

    gm_wgs84 = Constant(
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

    mass_wgs84 = Constant(
        abbrev='mass_wgs84',
        name='Mass of Earth for the WGS84 geodetic reference system, '
        'including the atmosphere',
        value=gm_wgs84.value / G.value,
        unit='kg',
        uncertainty=_np.sqrt((gm_wgs84.uncertainty / G.value)**2 +
                             (gm_wgs84.value * G.uncertainty / G.value**2)**2
                             ),
        reference='Derived from gm_wgs84 and G.')

    omega_wgs84 = Constant(
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

    gma_wgs84 = Constant(
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

    b_wgs84 = Constant(
        abbrev='b_wgs84',
        name='Semiminor axis of the WGS84 ellipsoid',
        value=6356752.3142,
        unit='m',
        uncertainty=0.0,
        reference='National Imagery and Mapping Agency (2000). Department of '
        'Defense World Geodetic System 1984: Its Definition and Relationship '
        'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
        'Mapping Agency.')

    r3_wgs84 = Constant(
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

    u0_wgs84 = Constant(
        abbrev='u0_wgs84',
        name="Theoretical normal gravity potential of the WGS84 ellipsoid",
        value=62636851.7146,
        unit='m2 / s2',
        uncertainty=0.0,
        reference='National Imagery and Mapping Agency (2000). Department of '
        'Defense World Geodetic System 1984: Its Definition and Relationship '
        'with Local Geodetic Systems. NIMA TR8350.2, National Imagery and '
        'Mapping Agency.')

    # == Moon ==

    gm_moon = Constant(
        abbrev='gm_moon',
        name='Gravitational constant times the mass of the Moon',
        value=4902.80007e9,
        unit='m3 / s2',
        uncertainty=0.00014e9,
        reference='Williams, J. G., A. S. Konopliv, D. H. Boggs, '
        'R. S. Park, D.-N. Yuan, F. G. Lemoine, S. Goossens, E. Mazarico, '
        'F. Nimmo, R. C. Weber, S. W. Asmar, H. J. Melosh, G. A. Neumann, '
        'R. J. Phillips, D. E. Smith, S. C. Solomon, M. M. Watkins, M. A. '
        'Wieczorek, J. C. Andrews-Hanna, J. W. Head, W. S. Kiefer, I. '
        'Matsuyama, P. J. McGovern, G. J. Taylor, and M. T. Zuber (2014). '
        'Lunar interior properties from the GRAIL mission, J. Geophys. Res. '
        'Planets, 119, 1546–1578, doi:10.1002/2013JE004559.')

    mass_moon = Constant(
        abbrev='mass_moon',
        name='Mass of the Moon',
        value=gm_moon.value / G.value,
        unit='kg',
        uncertainty=_np.sqrt((gm_moon.uncertainty / G.value)**2 +
                             (gm_moon.value * G.uncertainty / G.value**2)**2
                             ),
        reference='Derived from gm_moon and G.')

    r_moon = Constant(
        abbrev='r_moon',
        name='Mean radius of the Moon',
        value=1737151.0,
        unit='m',
        uncertainty=0.0,
        reference='LOLA2600p: Wieczorek, M. A. (2015). Gravity and topography '
        'of the terrestrial planets. In T. Spohn & G. Schubert (Eds.), '
        'Treatise on Geophysics, 2nd ed., Vol. 10, pp. 153–193). Oxford, '
        'Elsevier-Pergamon, doi:10.1016/B978-0-444-53802-4.00169-X.')

    density_moon = Constant(
        abbrev='density_moon',
        name='Mean density of the Moon',
        value=3 * mass_moon.value / (_np.pi * 4 * r_moon.value**3),
        unit='kg / m3',
        uncertainty=_np.sqrt((3 * mass_moon.uncertainty /
                             (_np.pi * 4 * r_moon.value**3))**2
                             + (3 * 3 * mass_moon.value *
                             r_moon.uncertainty /
                             (_np.pi * 4 * r_moon.value**4))**2
                             ),
        reference='Derived from mass_moon and r_moon.')

    g0_moon = Constant(
        abbrev='g0_moon',
        name='Surface gravity of the Moon, ignoring rotation and tides',
        value=gm_moon.value / r_moon.value**2,
        unit='m / s2',
        uncertainty=_np.sqrt((gm_moon.uncertainty / r_moon.value**2)**2
                             + (2 * gm_moon.value * r_moon.uncertainty
                             / r_moon.value**3)**2
                             ),
        reference='Derived from gm_moon and r_moon.')

    a_orbit_moon = Constant(
        abbrev='a_orbit_moon',
        name='Semimajor axis of the orbit of the Moon',
        value=384399.014e3,
        unit='m',
        uncertainty=0.0,
        reference='Williams, J. G., D. H. Boggs, and W. M. Folkner '
        '(2013), DE430 lunar orbit, physical librations, and surface '
        'coordinates, IOM 335-JW,DB, WF-20130722-016, July 22, 2013, '
        'Jet Propul. Lab., Pasadena, Calif.')

    omega_moon = Constant(
        abbrev='omega_moon',
        name='Angular spin rate of the Moon',
        value=2 * _np.pi / (27.321582 * 24 * 60 * 60),
        unit='rad / s',
        uncertainty=0.0,
        reference='Yoder, C. F. (1995). Astrometric and geodetic properties '
        'of Earth and the solar system. In: Ahrens TJ (ed.) Global Earth '
        'Physics: A Handbook of Physical Constants. AGU Reference Shelf, '
        'vol. 1, pp. 1-31. American Geophysical Union.')

    i_solid_moon = Constant(
        abbrev='i_solid_moon',
        name='Mean normalized moment of inertia of the solid portion of the '
        'Moon using the mean radius',
        value=0.393112,
        unit='',
        uncertainty=0.000012,
        reference='Williams, J. G., A. S. Konopliv, D. H. Boggs, '
        'R. S. Park, D.-N. Yuan, F. G. Lemoine, S. Goossens, E. Mazarico, '
        'F. Nimmo, R. C. Weber, S. W. Asmar, H. J. Melosh, G. A. Neumann, '
        'R. J. Phillips, D. E. Smith, S. C. Solomon, M. M. Watkins, M. A. '
        'Wieczorek, J. C. Andrews-Hanna, J. W. Head, W. S. Kiefer, I. '
        'Matsuyama, P. J. McGovern, G. J. Taylor, and M. T. Zuber (2014). '
        'Lunar interior properties from the GRAIL mission, J. Geophys. Res. '
        'Planets, 119, 1546–1578, doi:10.1002/2013JE004559.')

    beta_moon = Constant(
        abbrev='beta_moon',
        name='Libration parameter (C-A)/B of the Moon',
        value=631.0213e-6,
        unit='',
        uncertainty=0.0031e-6,
        reference='Williams, J. G., A. S. Konopliv, D. H. Boggs, '
        'R. S. Park, D.-N. Yuan, F. G. Lemoine, S. Goossens, E. Mazarico, '
        'F. Nimmo, R. C. Weber, S. W. Asmar, H. J. Melosh, G. A. Neumann, '
        'R. J. Phillips, D. E. Smith, S. C. Solomon, M. M. Watkins, M. A. '
        'Wieczorek, J. C. Andrews-Hanna, J. W. Head, W. S. Kiefer, I. '
        'Matsuyama, P. J. McGovern, G. J. Taylor, and M. T. Zuber (2014). '
        'Lunar interior properties from the GRAIL mission, J. Geophys. Res. '
        'Planets, 119, 1546–1578, doi:10.1002/2013JE004559.')

    gamma_moon = Constant(
        abbrev='beta_moon',
        name='Libration parameter (B-A)/C of the Moon',
        value=227.7317e-6,
        unit='',
        uncertainty=0.0042e-6,
        reference='Williams, J. G., A. S. Konopliv, D. H. Boggs, '
        'R. S. Park, D.-N. Yuan, F. G. Lemoine, S. Goossens, E. Mazarico, '
        'F. Nimmo, R. C. Weber, S. W. Asmar, H. J. Melosh, G. A. Neumann, '
        'R. J. Phillips, D. E. Smith, S. C. Solomon, M. M. Watkins, M. A. '
        'Wieczorek, J. C. Andrews-Hanna, J. W. Head, W. S. Kiefer, I. '
        'Matsuyama, P. J. McGovern, G. J. Taylor, and M. T. Zuber (2014). '
        'Lunar interior properties from the GRAIL mission, J. Geophys. Res. '
        'Planets, 119, 1546–1578, doi:10.1002/2013JE004559.')

    # == Mars ==

    gm_mars = Constant(
        abbrev='gm_mars',
        name='Gravitational constant times the mass of Mars',
        value=0.4282837581575610e+14,
        unit='m3 / s2',
        uncertainty=0.18167460e+6,
        reference='Konopliv A. S., R. S. Park ,W. M. Folkner (2016). '
        'An improved JPL Mars gravity field and orientation from Mars orbiter '
        'and lander tracking data, Icarus, 274, 253–260, '
        'doi:10.1016/j.icarus.2016.02.052')

    mass_mars = Constant(
        abbrev='mass_mars',
        name='Mass of Mars',
        value=gm_mars.value / G.value,
        unit='kg',
        uncertainty=_np.sqrt((gm_mars.uncertainty / G.value)**2 +
                             (gm_mars.value * G.uncertainty / G.value**2)**2
                             ),
        reference='Derived from gm_mars and G.')

    r_mars = Constant(
        abbrev='r_mars',
        name='Mean radius of Mars',
        value=3389.500e3,
        unit='m',
        uncertainty=0.0,
        reference='MarsTopo2600: Wieczorek, M. A. (2015). Gravity and '
        'topography of the terrestrial planets. In T. Spohn & G. Schubert '
        '(Eds.), Treatise on Geophysics, 2nd ed., Vol. 10, pp. 153–193). '
        'Oxford, Elsevier-Pergamon, doi:10.1016/B978-0-444-53802-4.00169-X.')

    density_mars = Constant(
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

    g0_mars = Constant(
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

    omega_mars = Constant(
        abbrev='omega_mars',
        name='Angular spin rate of Mars',
        value=350.891985307 * 2 * _np.pi / 360 / (24 * 60 * 60),
        unit='rad / s',
        uncertainty=0.000000003 * 2 * _np.pi / 360 / (24 * 60 * 60),
        reference='Konopliv A. S., R. S. Park ,W. M. Folkner (2016). '
        'An improved JPL Mars gravity field and orientation from Mars orbiter '
        'and lander tracking data, Icarus, 274, 253–260, '
        'doi:10.1016/j.icarus.2016.02.052')

    a_mars = Constant(
        abbrev='a_mars',
        name='Semimajor axis of the Mars reference ellipsoid',
        value=3395428.0,
        unit='m',
        uncertainty=19.0,
        reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
        'A new reference equipotential surface, and reference ellipsoid for '
        'the planet Mars. Earth, Moon, and Planets, 106, 1–13, '
        'doi:10.1007/s11038-009-9342-7.')

    b_mars = Constant(
        abbrev='b_mars',
        name='Semiminor axis of the Mars reference ellipsoid',
        value=3377678.0,
        unit='m',
        uncertainty=19.0,
        reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
        'A new reference equipotential surface, and reference ellipsoid for '
        'the planet Mars. Earth, Moon, and Planets, 106, 1–13, '
        'doi:10.1007/s11038-009-9342-7.')

    f_mars = Constant(
        abbrev='f_mars',
        name='Flattening of the Mars reference ellipsoid',
        value=(a_mars.value - b_mars.value) / a_mars.value,
        unit='',
        uncertainty=_np.sqrt((a_mars.uncertainty *
                             (a_mars.value - b_mars.value)
                             / a_mars.value**2)**2
                             + (_np.sqrt(a_mars.uncertainty**2 +
                             b_mars.uncertainty**2) / a_mars.value)**2
                             ),
        reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
        'A new reference equipotential surface, and reference ellipsoid for '
        'the planet Mars. Earth, Moon, and Planets, 106, 1–13, '
        'doi:10.1007/s11038-009-9342-7.')

    u0_mars = Constant(
        abbrev='u0_mars',
        name='Theoretical normal gravity potential of the reference ellipsoid',
        value=12654875.0,
        unit='m2 / s2',
        uncertainty=69.0,
        reference='Ardalan A. A, R. Karimi, and E. W. Grafarend (2010). '
        'A new reference equipotential surface, and reference ellipsoid for '
        'the planet Mars. Earth, Moon, and Planets, 106, 1–13, '
        'doi:10.1007/s11038-009-9342-7.')

except ImportError:
    raise ImportError('To use the pyshtools constant subpackage, you must '
                      'install astropy.')
