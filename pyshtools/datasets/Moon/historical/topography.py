'''
Historical datasets related to the shape of Earth's Moon.

Topography
----------
GLTM2B (Clementine)       :  Smith et al. (1997)
MoonTopo2600p (LRO)       :  Wieczorek (2015)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import DOIDownloader as _DOIDownloader
from ....shclasses import SHCoeffs as _SHCoeffs


def MoonTopo2600p(lmax=2600):
    '''
    MoonTopo2600p is a 2600 degree and order spherical harmonic model of the
    shape of Earth's Moon in a principal axis coordinate system. The
    coefficients are in units of meters.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Wieczorek, M.A. (2015). Gravity and Topography of the Terrestrial Planets,
        Treatise on Geophysics, 2nd edition, Oxford, 153-193,
        doi:10.1016/B978-0-444-53802-4.00169-X.
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.3870924/MoonTopo2600p.shape.gz",
        known_hash="sha256:193146df894e2fef796df9d6142c78fae6fa5c183fd79d3f79eeb356602af69a",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='MoonTopo2600p',
                               units='m', encoding='utf-8')


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


__all__ = ['GLTM2B', 'MoonTopo2600p']
