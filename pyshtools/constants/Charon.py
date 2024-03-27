"""
pyshtools constants for Pluto's moon Charon.

Each object is an astropy Constant that possesses the attributes name, value,
unit, uncertainty, and reference.
"""
import numpy as _np

from astropy.constants import Constant as _Constant
from astropy.constants import G as _G


gm = _Constant(
    abbrev='gm_charon',
    name='Gravitational constant times the mass of Charon',
    value=105.88e9,
    unit='m3 / s2',
    uncertainty=1.0e9,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

mass = _Constant(
    abbrev='mass_charon',
    name='Mass of Charon',
    value=gm.value / _G.value,
    unit='kg',
    uncertainty=_np.sqrt((gm.uncertainty / _G.value)**2 +
                         (gm.value * _G.uncertainty / _G.value**2)**2
                         ),
    reference='Derived from gm_charon and G.')

mean_radius = _Constant(
    abbrev='mean_radius_charon',
    name='Mean radius of Charon',
    value=606.0e3,
    unit='m',
    uncertainty=1.0e3,
    reference='Nimmo, F., Umurhan, O., Lisse, C. M., Bierson, C. J., Lauer, '
    'T. R., Buie, M. W., Throop, H. B., Kammer, J. A., Roberts, J. H., '
    'McKinnon, W. B., Zangari, A. M., Moore, J. M., Stern, S. A., Young, '
    'L. A., Weaver, H. A., Olkin, C. B., & Ennico, K. (2017). Mean radius and '
    'shape of Pluto and Charon from New Horizons images. Icarus, 287, 12–29. '
    'https://doi.org/10.1016/j.icarus.2016.06.027')

r = mean_radius

volume_equivalent_radius = _Constant(
    abbrev='volume_equivalent_radius_charon',
    name='Volume equivalent radius of Charon',
    value=mean_radius.value,
    unit='m',
    uncertainty=mean_radius.uncertainty,
    reference='Equal to mean_radius_charon')

volume = _Constant(
    abbrev='volume_charon',
    name='Volume of Charon',
    value=(4 * _np.pi / 3) * volume_equivalent_radius.value**3,
    unit='m',
    uncertainty=(8 * _np.pi / 3) * volume_equivalent_radius.value**2 *
    volume_equivalent_radius.uncertainty,
    reference='Derived from volume_equivalent_radius_charon')

mean_density = _Constant(
    abbrev='mean_density_charon',
    name='Mean density of Charon',
    value=3 * mass.value / (_np.pi * 4 * volume_equivalent_radius.value**3),
    unit='kg / m3',
    uncertainty=_np.sqrt((3 * mass.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**3))**2
                         + (3 * 3 * mass.value *
                         volume_equivalent_radius.uncertainty /
                         (_np.pi * 4 * volume_equivalent_radius.value**4))**2
                         ),
    reference='Derived from mass_charon and volume_equivalent_radius_charon.')

gravity_mean_radius = _Constant(
    abbrev='gravity_mean_radius_charon',
    name='Surface gravity of Charon, ignoring rotation and tides',
    value=gm.value / mean_radius.value**2,
    unit='m / s2',
    uncertainty=_np.sqrt((gm.uncertainty / mean_radius.value**2)**2
                         + (2 * gm.value * mean_radius.uncertainty
                         / mean_radius.value**3)**2
                         ),
    reference='Derived from gm_charon and mean_radius_charon.')

angular_velocity = _Constant(
    abbrev='angular_velocity_charon',
    name='Angular spin rate of Charon',
    value=1.1385591834674098e-05,
    unit='rad / s',
    uncertainty=0.,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

rotational_period = _Constant(
    abbrev='rotational_period_charon',
    name='Rotational period of Charon',
    value=2. * _np.pi / angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from angular_velocity_charon')

orbit_semimajor_axis = _Constant(
    abbrev='orbit_semimajor_axis_charon',
    name='Semimajor axis of the orbit of Charon about Pluto',
    value=19596.e3,
    unit='m',
    uncertainty=0.,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

orbit_eccentricity = _Constant(
    abbrev='orbit_eccentricity_charon',
    name='Eccentricity of the orbit of Charon about Pluto',
    value=0.00005,
    unit='',
    uncertainty=0.,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

orbit_inclination = _Constant(
    abbrev='orbit_inclination_charon',
    name='Inclination of the orbit of Charon about Pluto',
    value=0.0,
    unit='degrees',
    uncertainty=0.,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

orbit_angular_velocity = _Constant(
    abbrev='orbit_angular_velocity_charon',
    name='Orbital angular velocity of Charon',
    value=1.1385591834674098e-05,
    unit='rad / s',
    uncertainty=0.,
    reference='Brozović, M., Showalter, M. R., Jacobson, R. A., & Buie, M. W. '
    '(2015). The orbits and masses of satellites of Pluto. Icarus, 246, '
    '317–329. https://doi.org/10.1016/j.icarus.2014.03.015')

orbit_period = _Constant(
    abbrev='orbit_period_charon',
    name='Orbital period of Charon',
    value=2. * _np.pi / orbit_angular_velocity.value,
    unit='s',
    uncertainty=2. * _np.pi * orbit_angular_velocity.uncertainty /
    angular_velocity.value**2,
    reference='Derived from orbit_angular_velocity_charon')

__all__ = ['gm', 'mass', 'mean_radius', 'r', 'volume_equivalent_radius',
           'volume', 'gravity_mean_radius', 'mean_density', 'angular_velocity',
           'rotational_period', 'orbit_semimajor_axis', 'orbit_eccentricity',
           'orbit_inclination', 'orbit_angular_velocity', 'orbit_period']
