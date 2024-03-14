'''
Datasets related to Jupiter's moon Io.

Gravity
-------
Anderson2001  :  Anderson et al. (2001)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Io import omega as _omega


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
    Anderson, J. D., Jacobson, R. A., Lau, E. L., Moore, W. B., & Schubert, G.
        (2001). Io's gravity field and interior structure. J. Geophys. Res.,
        106, 32963–32969. https://doi.org/10.1029/2000JE001367
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10817282/Anderson2001_Io_gravity.sh",  # noqa: E501
        known_hash="sha256:d2090309dc67ffa2c2c2e96eb9a3cb8a86b3631be39f43bd4112e07367cdfeaa",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Anderson2001 (Io)',
                                   encoding='utf-8', omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Anderson2001']
