"""
pyshtools constants for Jupiter's moon Io.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_europa',
    name='Gravitational constant times the mass of europa',
    value=3202.72e9,
    unit='m3 / s2',
    uncertainty=0.05e9,
    reference='Anderson, J. D., Schubert, G., Jacobson, R. A., Lau, E. L., '
    "Moore, W. B., & Sjogren, W. L. (1998). Europa's differentiated internal "
    'structure: Inferences from four Galileo encounters. Science, 281, '
    '2019–2022. https://doi.org/10.1126/science.281.5385.2019')

mass = _Constant(
    abbrev='mass_europa',
    name='Mass of Europa',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_Europa and G.')

mean_radius = _Constant(
    abbrev='mean_radius_europa',
    name='Mean radius of Europa',
    value=1560.8e3,
    unit='m',
    uncertainty=0.3e3,
    reference='Nimmo, F., Thomas, P., Pappalardo, R., & Moore, W. (2007). The '
    'global shape of Europa: Constraints on lateral shell thickness '
    'variations. Icarus, 191(1), 183–192. '
    'https://doi.org/10.1016/j.icarus.2007.04.021')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_europa',
    name='Volume equivalent radius of Europa',
    value=mean_radius.value,
    unit='m',
    uncertainty=mean_radius.uncertainty,
    reference='Equal to mean_radius_europa')

volume = _Constant(
    abbrev='volume_europa',
    name='Volume of Europa',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_europa')

mean_density = _Constant(
    abbrev='mean_density_europa',
    name='Mean density of Europa',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_europa and volume_equivalent_radius_europa.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_europa',
    name='Gravity at mean radius of Europa, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_europa and mean_radius_europa.')

omega = _Constant(
    abbrev='omega_europa',
    name='Angular spin rate of Europa',
    value=101.3747235 * 2. * _np.pi / 360. / (24. * 60. * 60.),
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
