'''
Historical datasets related to the shape of Mars.

Shape
-----
MarsTopo2600     :  Wieczorek (2015)
MarsTopo719      :  Wieczorek (2015)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ....shclasses import SHCoeffs as _SHCoeffs


def MarsTopo2600(lmax=2600):
    '''
    MarsTopo2600 is a 2600 degree and order spherical harmonic model of the
    shape of the planet Mars. The coefficients are in units of meters.

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
        url="doi:10.5281/zenodo.3870922/MarsTopo2600.shape.gz",
        known_hash="sha256:8882a9ee7ee405d971b752028409f69bd934ba5411f1c64eaacd149e3b8642af",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='MarsTopo2600',
                               units='m', encoding='utf-8')


def MarsTopo719(lmax=719):
    '''
    MarsTopo719 is a 719 degree and order spherical harmonic model of the
    shape of the planet Mars. The coefficients are in units of meters. This
    dataset is a truncated version of MarsTopo2600.

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
        url="doi:10.5281/zenodo.6475460/MarsTopo719.shape.gz",
        known_hash="sha256:37a98efae5eab7c85260f4b43315fe9fcf44247a61581bed1b6f7f10f79adea0",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='MarsTopo719',
                               units='m', encoding='utf-8')


__all__ = ['MarsTopo2600', 'MarsTopo719']
