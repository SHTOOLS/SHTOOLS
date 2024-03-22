"""
pyshtools constants for the asteroid (4) Vesta.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G

gm = _Constant(
    abbrev='gm_vesta',
    name='Gravitational constant times the mass of (4) Vesta',
    value=317288244969.3,
    unit='m3 / s2',
    uncertainty=4064.9,
    reference='VESTA20H: Konopliv, A.S., Asmar, S.W., Park, R.S., Bills, '
    'B.G., Centinello, F., Chamberlin, A.B., Ermakov, A., Gaskell, R.W., '
    'Rambaux, N., Raymond, C.A., Russell, C.T., Smith, D.E., Tricarico, P., '
    'Zuber, M.T. (2014). The Vesta gravity field, spin pole and rotation '
    'period, landmark positions, and ephemeris from the Dawn tracking and '
    'optical data. Icarus, 240, 103-117, doi:10.1016/j.icarus.2013.09.005.')

mass = _Constant(
    abbrev='mass_vesta',
    name='Mass of Vesta',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_vesta and G.')

mean_radius = _Constant(
    abbrev='mean_radius_vesta',
    name='Mean radius of (4) Vesta',
    value=260362.3,
    unit='m',
    uncertainty=0.0,
    reference='DLR_SPG_shape: Wieczorek, M. (2024). Spherical harmonic models '
    'of the shape of the asteroid (4) Vesta [DLR SPG] (1.0.0) [Data set]. '
    'Zenodo. https://doi.org/10.5281/zenodo.10800929')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_vesta',
    name='Volume equivalent radius of (4) Vesta',
    value=261590.3,
    unit='m',
    uncertainty=0.,
    reference='Computed using DLR_SPG_shape and SHCoeffs.volume()')

volume = _Constant(
    abbrev='volume_vesta',
    name='Volume of (4) Vesta',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_vesta')

mean_density = _Constant(
    abbrev='mean_density_vesta',
    name='Mean density of (4) Vesta',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_vesta and volume_equivalent_radius_vesta.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_vesta',
    name='Gravity at the mean radius of (4) Vesta, ignoring rotation',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_vesta and mean_radius_vesta.')

angular_velocity = _Constant(
    abbrev='angular_velocity_vesta',
    name='Angular spin rate of (4) Vesta',
    value=1617.3331279 * 2. * _np.pi / 360. / (24. * 60. * 60.),
    unit='rad / s',
    uncertainty=0.,
    reference='VESTA20H: Konopliv, A.S., Asmar, S.W., Park, R.S., Bills, '
    'B.G., Centinello, F., Chamberlin, A.B., Ermakov, A., Gaskell, R.W., '
    'Rambaux, N., Raymond, C.A., Russell, C.T., Smith, D.E., Tricarico, P., '
    'Zuber, M.T. (2014). The Vesta gravity field, spin pole and rotation '
    'period, landmark positions, and ephemeris from the Dawn tracking and '
    'optical data. Icarus, 240, 103-117, doi:10.1016/j.icarus.2013.09.005.')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity']
