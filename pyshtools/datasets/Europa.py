'''
Datasets related to Jupiter's moon Europa.

Gravity
-------
Anderson1998  :  Anderson et al. (1998)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Europa import omega as _omega


def Anderson1998(lmax=2):
    '''
    Anderson1998 is a JPL spherical harmonic model of the gravitational
    potential of Europa to degree and order 2.

    Parameters
    ----------
    lmax : int, optional, default = 2
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Anderson, J. D., Schubert, G., Jacobson, R. A., Lau, E. L., Moore, W. B.,
        & Sjogren, W. L. (1998). Europa's differentiated internal structure:
        Inferences from four Galileo encounters. Science, 281, 2019–2022.
        https://doi.org/10.1126/science.281.5385.2019
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10817282/Anderson1998_Europa_gravity.sh",  # noqa: E501
        known_hash="sha256:52d5f62ed31fd1dce8324fdcee8eebe105509bbeafaa254a9056a0eb46314615",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Anderson1998 (Europa)',
                                   encoding='utf-8', omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Anderson1998']
