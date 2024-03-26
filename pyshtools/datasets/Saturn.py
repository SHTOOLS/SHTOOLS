'''
Datasets related to the planet Saturn.

Gravity
-------
Jacobson2022_gravity  :  Jacobson et al. (2022)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import create as _create
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHCoeffs as _SHCoeffs
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Saturn import angular_velocity as _omega


def Jacobson2022_gravity(lmax=12):
    '''
    Jacobson2022_gravity is a degree and order 12 spherical harmonic model of the
    gravitational potential of Saturn.

    Parameters
    ----------
    lmax : int, optional, default = 3
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Jacobson, R. (2022). The Orbits of the Main Saturnian Satellites, the
        Saturnian System Gravity Field, and the Orientation of Saturn's Pole.
        The Astronomical Journal, 164, 199.
        https://doi.org/10.3847/1538-3881/ac90c9
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10878627/Jacobson2022.sh",
        known_hash="sha256:9f39f74cca125e3c1e44a9fded2ddb648397b1c13850b2b136465c18589174e6",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Jacobson2022_gravity (Saturn)',
                                   encoding='utf-8',
                                   omega=_omega.value,
                                   normalization='unnorm')

__all__ = ['Jacobson2022_gravity']
