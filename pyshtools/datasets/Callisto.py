'''
Datasets related to Jupiter's moon Callisto.

Gravity
-------
Anderson2001  :  Anderson et al. (2001)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Callisto import angular_velocity as _omega


def Anderson2001(lmax=2):
    '''
    Anderson2001 is a JPL spherical harmonic model of the gravitational
    potential of Io to degree and order 2.

    Parameters
    ----------
    lmax : int, optional, default = 2
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Anderson, J. D., Jacobson, R. A., McElrath, T. P., Moore, W. B., &
        Schubert, G. (2001). Shape, mean radius, gravity field, and interior
        structure of Callisto. Icarus, 153(1), 157–161.
        https://doi.org/10.1006/icar.2001.6664
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10817282/Anderson2001_Callisto_gravity.sh",
        known_hash="sha256:472d3c6837e96c2aa89c2fb154b71b187248c01e0cb47482d0e743cd02002c6b",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Anderson2001 (Callisto)',
                                   encoding='utf-8', omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Anderson2001']
