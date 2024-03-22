"""
pyshtools constants for Earth's Moon.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
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
    'Planets, 119, 1546-1578, doi:10.1002/2013JE004559.')

mass = _Constant(
    abbrev='mass_moon',
    name='Mass of the Moon',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_moon and G.')

mean_radius = _Constant(
    abbrev='mean_radius_moon',
    name='Mean radius of the Moon',
    value=1737151.,
    unit='m',
    uncertainty=0.,
    reference='LOLA_shape: Wieczorek, M. A. (2024). Spherical harmonic models'
    'of the shape of the Moon (1.0.0) [Data set]. Zenodo. '
    'https://doi.org/10.5281/zenodo.10796823')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_moon',
    name='Volume equivalent radius of the Moon',
    value=1737154.4,
    unit='m',
    uncertainty=0.,
    reference='Computed using LOLA_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_moon',
    name='Volume of Moon',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_moon')

mean_density = _Constant(
    abbrev='mean_density_moon',
    name='Mean density of the Moon',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_moon and volume_equivalent_radius_moon.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_moon',
    name='Gravity at the mean radius of the Moon, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_moon and mean_radius_moon.')

a_orbit = _Constant(
    abbrev='a_orbit_moon',
    name='Semimajor axis of the orbit of the Moon',
    value=384399.014e3,
    unit='m',
    uncertainty=0.0,
    reference='Williams, J. G., D. H. Boggs, and W. M. Folkner '
    '(2013), DE430 lunar orbit, physical librations, and surface '
    'coordinates, IOM 335-JW,DB, WF-20130722-016, July 22, 2013, '
    'Jet Propul. Lab., Pasadena, Calif.')

angular_velocity = _Constant(
    abbrev='angular_velocity_moon',
    name='Angular spin rate of the Moon',
    value=2 * _np.pi / (27.321582 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.0,
    reference='Yoder, C. F. (1995). Astrometric and geodetic properties '
    'of Earth and the solar system. In: Ahrens TJ (ed.) Global Earth '
    'Physics: A Handbook of Physical Constants. AGU Reference Shelf, '
    'vol. 1, pp. 1-31. American Geophysical Union.')

i_solid = _Constant(
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
    'Planets, 119, 1546-1578, doi:10.1002/2013JE004559.')

beta = _Constant(
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
    'Planets, 119, 1546-1578, doi:10.1002/2013JE004559.')

gamma = _Constant(
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
    'Planets, 119, 1546-1578, doi:10.1002/2013JE004559.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'mean_density', 'gravity_mean_radius', 'a_orbit',
           'angular_velocity', 'i_solid', 'beta', 'gamma']
