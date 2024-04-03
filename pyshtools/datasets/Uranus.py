'''
Datasets related to the planet Uranus.

Gravity
-------
Jacobson2014_gravity  :  Jacobson et al. (2014)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Uranus import angular_velocity as _omega


def Jacobson2014_gravity(lmax=6):
    '''
    Jacobson2014_gravity is a degree 6 and order 0 spherical harmonic model of
    the gravitational potential of Uranus.

    Parameters
    ----------
    lmax : int, optional, default = 6
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Jacobson, R. A. (2014). The orbits of the Uranian satellites and rings,
        the gravity field of the Uranian system, and the orientation of the
        pole of Uranus. The Astronomical Journal, 148(5), 76.
        https://doi.org/10.1088/0004-6256/148/5/76
    '''
    if lmax < 0:
        lmax = 6

    fname = _retrieve(
        url="doi:10.5281/zenodo.10883631/Jacobson2014.sh",
        known_hash="sha256:18598990169d75375139d6e195ccdf4b9b784a444c5a8606bc495b10495541b9",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Jacobson2014_gravity (Uranus)',
                                   encoding='utf-8',
                                   omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Jacobson2014_gravity']
