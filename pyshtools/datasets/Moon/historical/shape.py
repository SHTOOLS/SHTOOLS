'''
Historical datasets related to the shape of Earth's Moon.

Shape
-----
GLTM2B (Clementine)       :  Smith et al. (1997)
MoonTopo2600p (LRO)       :  Wieczorek (2015)
LOLA_shape_pa (LRO)       :  Wieczorek (2024)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import DOIDownloader as _DOIDownloader
from ....shclasses import SHCoeffs as _SHCoeffs


def LOLA_shape_pa(lmax=719):
    '''
    LOLA_shape_pa is a spherical harmonic model of the shape of Earth's Moon in
    a principal axis coordinate system based on LOLA laser altimetry data
    obtained by the Lunar Reconaissance Orbiter mission. The maximum spherical
    harmonic degree of the model is 5759, which has an effective spatial
    resolution of 64 pixels per degree. Three lower resolution models are
    available in this archive (with lmax of 719, 1439 and 2879), and only the
    smallest that is required by the user input lmax will be downloaded. If
    lmax is not specified, the lowest resolution model (719) will be returned.
    If a negative value for lmax is specified, the maximum resolution model
    will be returned. The coefficients are in units of meters.

    This shape model uses the same coordinate system as most lunar gravity
    models. For a mean Earth/polar axis model, use LOLA_shape instead.

    Parameters
    ----------
    lmax : int, optional, default = 719
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Wieczorek, M. (2024). Spherical harmonic models of the shape of the Moon
        (principal axis coordinate system) [LOLA] (1.0.1) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.10820750
    LRO LOLA Team (2013). LRO-L-LOLA-4-GDR-V1.0, NASA Planetary Data System.
    '''
    archive = _create(
        path=_os_cache('pyshtools'),
        base_url="doi:10.5281/zenodo.10820750",
        registry={
            "Moon_LOLA_shape_pa_5759.bshc.gz": "sha256:1569338a88475e184ac8b7424327b0e03cb1f8f0a3cec70bd9bfe41635f68671",  # noqa: E501
            "Moon_LOLA_shape_pa_2879.bshc.gz": "sha256:6ee87880956fdbdf85f9c0b2e9996514735198fbcb19fd0d0c7881bf50c50496",  # noqa: E501
            "Moon_LOLA_shape_pa_1439.bshc.gz": "sha256:6c4619f845c9902d999879cbf6956c368290a17fc2887d41267800aade386c56",  # noqa: E501
            "Moon_LOLA_shape_pa_719.bshc.gz": "sha256:71877e8c1dd80205941b6ca0e7df73943abc52be125d1d8cc76bea2dcee5942b",  # noqa: E501
            },
        )

    if lmax < 0:
        lmax = 5759

    if lmax >= 0 and lmax <= 719:
        fname = archive.fetch("Moon_LOLA_shape_pa_719.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 719 and lmax <= 1439:
        fname = archive.fetch("Moon_LOLA_shape_pa_1439.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 1439 and lmax <= 2879:
        fname = archive.fetch("Moon_LOLA_shape_pa_2879.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
    else:
        fname = archive.fetch("Moon_LOLA_shape_pa_5759.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
        lmax = min(lmax, 5759)

    return _SHCoeffs.from_file(fname, lmax=lmax, name='LOLA_shape_pa (Moon)',
                               units='m', format='bshc')


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
    if lmax < 0:
        lmax = 2600

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
    if lmax < 0:
        lmax = 70

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


__all__ = ['GLTM2B', 'MoonTopo2600p', 'LOLA_shape_pa']
