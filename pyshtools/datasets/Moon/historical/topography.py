'''
Historical datasets related to the shape of Earth's Moon.

Topography
-------
GLTM2B (Clementine)       :  Smith et al. (1997)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ....shclasses import SHCoeffs as _SHCoeffs


def GLTM2B(lmax=70):
    '''
    GLTM-2B is a GSFC 70 degree and order spherical harmonic model of the
    shape of the Moon derived from Clementine laser altimetry data.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Smith, D. E., Zuber, M. T., Neumann, G. A., Lemoine, F. G. (1997).
        Topography of the Moon from the Clementine lidar. Journal of
        Geophysical Research, 102 (E1), 1591-1611, doi:10.1029/96JE02940.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/lunar/clem1-gravity-topo-v1/cl_8xxx/topo/gltm2bsh.tab",  # noqa: E501
        known_hash="sha256:002bc2bccaaef6175f3e42f717e06dd267cd699d6d5bae6e6c23d344ca7293d4",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    clm = _SHCoeffs.from_file(fname, lmax=lmax, header=False, errors=False,
                              name='GLTM-2B', encoding='utf-8')
    clm += 1738.e3
    return clm


__all__ = ['GLTM2B']
