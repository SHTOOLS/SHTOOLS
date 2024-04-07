'''
Datasets related to the planet Jupiter.

Gravity
-------
Kaspi2023_gravity  :  Kaspi et al. (2023)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Jupiter import angular_velocity as _omega


def Kaspi2023_gravity(lmax=40):
    '''
    Kaspi2023_gravity is a degree 40 and order 0 spherical harmonic model of
    the gravitational potential of Jupiter.

    Parameters
    ----------
    lmax : int, optional, default = 40
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Kaspi, Y., Galanti, E., Park, R. S., Duer, K., Gavriel, N., Durante, D.,
        Iess, L., Parisi, M., Buccino, D. R., Guillot, T., Stevenson, D. J.,
        & Bolton, S. J. (2023). Observational evidence for cylindrically
        oriented zonal flows on Jupiter. Nature Astronomy, 7(12), 1463–1472.
        https://doi.org/10.1038/s41550-023-02077-8
    '''
    if lmax < 0:
        lmax = 40

    fname = _retrieve(
        url="doi:10.5281/zenodo.10886089/Kaspi2023.sh",
        known_hash="sha256:9517c67f623c5e252ce77ef2650f8efed798ac7b7fc68768615e610a3ac34db1",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='kKaspi2023_gravity (Jupiter)',
                                   encoding='utf-8',
                                   omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Kaspi2023_gravity']
