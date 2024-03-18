"""
pyshtools constants for the planet Venus.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_venus',
    name='Gravitational constant times the mass of Venus',
    value=324858592079000.,
    unit='m3 / s2',
    uncertainty=6376000.0,
    reference='MGNP180U: Konopliv A. S., W. B. Banerdt, and W. L. Sjogren '
    '(1999) Venus gravity: 180th degree and order model. Icarus 139: 3-18.'
    'doi:10.1006/icar.1999.6086.')

mass = _Constant(
    abbrev='mass_venus',
    name='Mass of Venus',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_venus and G.')

mean_radius = _Constant(
    abbrev='mean_radius_venus',
    name='Mean radius of Venus',
    value=6051877.4,
    unit='m',
    uncertainty=0.0,
    reference='VenusTopo719: Wieczorek, M. A. (2015). Gravity and '
    'topography of the terrestrial planets. In T. Spohn & G. Schubert '
    '(Eds.), Treatise on Geophysics, 2nd ed., Vol. 10, pp. 153-193). '
    'Oxford, Elsevier-Pergamon, doi:10.1016/B978-0-444-53802-4.00169-X.')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_venus',
    name='Volume equivalent radius of Venus',
    value=6051877.5,
    unit='m',
    uncertainty=0.,
    reference='Computed using VenusTopo719 and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_venus',
    name='Volume of Venus',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_venus')

mean_density = _Constant(
    abbrev='mean_density_venus',
    name='Mean density of Venus',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_venus and volume_equivalent_radius_venus.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_venus',
    name='Gravity at the mean radius of Venus, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_venus and mean_radius_venus.')

omega = _Constant(
    abbrev='omega_venus',
    name='Angular spin rate of Venus',
    value=-2 * _np.pi / (243.0226 * 24 * 60 * 60),
    unit='rad / s',
    uncertainty=0.0013 * 2 * _np.pi / (243.0226**2 * 24 * 60 * 60),
    reference='Margot, J.-L., D. B. Campbell, J. D. Giorgini, J. S. Jao, '
    'L. G. Snedeker, F. D. Ghigo and A. Bonsall, Spin state and moment of '
    'inertia of Venus (2021), Nature Astronomy, '
    'doi:10.1038/s41550-021-01339-7.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'omega']
