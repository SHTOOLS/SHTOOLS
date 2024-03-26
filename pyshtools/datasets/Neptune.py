'''
Datasets related to the planet Neptune.

Gravity
-------
Jacobson2009_gravity  :  Jacobson et al. (2009)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Neptune import angular_velocity as _omega


def Jacobson2009_gravity(lmax=4):
    '''
    Jacobson2009_gravity is a degree 4 and order 0 spherical harmonic model of
    the gravitational potential of Neptune.

    Parameters
    ----------
    lmax : int, optional, default = 4
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Jacobson, R. A. (2009). The orbits of the Neptunian satellites and the
        orientation of the pole of Neptune. The Astronomical Journal, 137(5),
        4322. https://doi.org/10.1088/0004-6256/137/5/4322

    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10883883/Jacobson2009.sh",
        known_hash="sha256:200d6e00446edfa5a25f82646a5d0e882d36c54dae862b20411d90b4d20ad5de",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Jacobson2009_gravity (Neptune)',
                                   encoding='utf-8',
                                   omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Jacobson2009_gravity']
