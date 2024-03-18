"""
pyshtools constants for Jupiter's moon Io.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_io',
    name='Gravitational constant times the mass of Io',
    value=5959.91e9,
    unit='m3 / s2',
    uncertainty=0.28e09,
    reference='Anderson, J. D., Jacobson, R. A., Lau, E. L., Moore, W. B., & '
    "Schubert, G. (2001). Io's gravity field and interior structure. J. "
    'Geophys. Res., 106, 32963–32969. https://doi.org/10.1029/2000JE001367')

mass = _Constant(
    abbrev='mass_io',
    name='Mass of Io',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_io and G.')

mean_radius = _Constant(
    abbrev='mean_radius_io',
    name='Mean radius of Io',
    value=1821.49e3,
    unit='m',
    uncertainty=500.,
    reference='Thomas, P. C., Davies, M. E., Colvin, T. R., Oberst, J., '
    'Schuster, P., Neukum, G., Carr, M. H., McEwen, A., Schubert, G., & '
    'Belton, M. J. S. (1998). The Shape of Io from Galileo Limb Measurements. '
    'Icarus, 135(1), 175–180. https://doi.org/10.1006/icar.1998.5987')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_io',
    name='Volume equivalent radius of io',
    value=mean_radius.value,
    unit='m',
    uncertainty=mean_radius.uncertainty,
    reference='Equal to mean_radius_io')

volume = _Constant(
    abbrev='volume_io',
    name='Volume of Io',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_io')

mean_density = _Constant(
    abbrev='mean_density_io',
    name='Mean density of Io',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_io and volume_equivalent_radius_io.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_io',
    name='Surface gravity of Io, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_io and mean_radius_io.')

omega = _Constant(
    abbrev='omega_Io',
    name='Angular spin rate of Io',
    value=203.4889538 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='Archinal, B. A., Acton, C. H., A’Hearn, M. F., Conrad, A., '
    'Consolmagno, G. J., Duxbury, T., Hestroffer, D., Hilton, J. L., Kirk, '
    'R. L., Klioner, S. A., McCarthy, D., Meech, K., Oberst, J., Ping, J., '
    'Seidelmann, P. K., Tholen, D. J., Thomas, P. C., & Williams, I. P. '
    '(2018). Report of the IAU Working Group on Cartographic Coordinates and '
    'Rotational Elements: 2015. Celestial Mechanics and Dynamical Astronomy, '
    '130(3), 22. https://doi.org/10.1007/s10569-017-9805-5')


__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'omega']
