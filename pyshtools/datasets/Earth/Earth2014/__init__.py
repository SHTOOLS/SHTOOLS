'''
Earth2014: Topography of Earth with respect to mean sea level, expanded to
degree and order 2160.

surface   :  Earth's physical surface (including water and ice masses)
bedrock   :  Earth's bedrock (excluding water and ice masses)
tbi       :  Topography of Earth's bedrock and ice masses (excluding water)
ret       :  Earth's rock-equivalent topography
ice       :  Thickness of Earth's ice sheets

Reference
---------
Hirt, C., Rexer, M. (2015). Earth2014: 1 arc-min shape, topography, bedrock
    and ice-sheet models - available as gridded data and degree-10,800
    spherical harmonics, International Journal of Applied Earth Observation and
    Geoinformation, 39, 103-112, doi:10.10.1016/j.jag.2015.03.001.
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ....shclasses import SHCoeffs as _SHCoeffs


def surface(lmax=2160):
    '''
    Earth's physical surface: Harmonic model of the interface between the
    atmosphere and Earth's surface (including water and ice).

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Rexer, M. (2015). Earth2014: 1 arc-min shape, topography, bedrock
        and ice-sheet models - available as gridded data and degree-10,800
        spherical harmonics, International Journal of Applied Earth Observation
        and Geoinformation, 39, 103-112, doi:10.10.1016/j.jag.2015.03.001.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2014/data_5min/shcs_to2160/Earth2014.SUR2014.degree2160.bshc",  # noqa: E501
        known_hash="sha256:5694260d135c17427270ed18d48af23f4788e5fbc1dfb9dcb19f1cf9b401c9ce",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, format='bshc', units='m')


def bedrock(lmax=2160):
    '''
    Earth's bedrock: Harmonic model of the surface relief, excluding water
    and ice masses.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Rexer, M. (2015). Earth2014: 1 arc-min shape, topography, bedrock
        and ice-sheet models - available as gridded data and degree-10,800
        spherical harmonics, International Journal of Applied Earth Observation
        and Geoinformation, 39, 103-112, doi:10.10.1016/j.jag.2015.03.001.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2014/data_5min/shcs_to2160/Earth2014.BED2014.degree2160.bshc",  # noqa: E501
        known_hash="sha256:146dcc80f17d201352d391aa90f487f5ed16006a6a3966add2d023b998727af7",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, format='bshc', units='m')


def tbi(lmax=2160):
    '''
    Topography of Earth's bedrock and ice: Harmonic model of the surface
    relief, excluding water masses.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Rexer, M. (2015). Earth2014: 1 arc-min shape, topography, bedrock
        and ice-sheet models - available as gridded data and degree-10,800
        spherical harmonics, International Journal of Applied Earth Observation
        and Geoinformation, 39, 103-112, doi:10.10.1016/j.jag.2015.03.001.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2014/data_5min/shcs_to2160/Earth2014.TBI2014.degree2160.bshc",  # noqa: E501
        known_hash="sha256:84a72ef25fe26fd746d2c6988c01840a12e9e414ee266e9357acb81faaaa6d5f",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, format='bshc', units='m')


def ret(lmax=2160):
    '''
    Earth's rock-equivalent topography: Harmonic model of Earth's rock
    equivalent topography, where ice and water masses are condensed to layers
    of rock using a density of 2670 kg/m3.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Rexer, M. (2015). Earth2014: 1 arc-min shape, topography, bedrock
        and ice-sheet models - available as gridded data and degree-10,800
        spherical harmonics, International Journal of Applied Earth Observation
        and Geoinformation, 39, 103-112, doi:10.10.1016/j.jag.2015.03.001.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2014/data_5min/shcs_to2160/Earth2014.RET2014.degree2160.bshc",  # noqa: E501
        known_hash="sha256:1fa87749532811614c00e9723dae1aa312e5511570d101bccb74bc40cb7dd5d1",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, format='bshc', units='m')


def ice(lmax=2160):
    '''
    Thickness of Earth's ice sheets: Harmonic model of the thickness of
    Earth's ice masses.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Rexer, M. (2015). Earth2014: 1 arc-min shape, topography, bedrock
        and ice-sheet models - available as gridded data and degree-10,800
        spherical harmonics, International Journal of Applied Earth Observation
        and Geoinformation, 39, 103-112, doi:10.10.1016/j.jag.2015.03.001.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2014/data_5min/shcs_to2160/Earth2014.ICE2014.degree2160.bshc",  # noqa: E501
        known_hash="sha256:04cd185cc668eba6f9bd1db527c4985703cce8d1a9fb993509e625e2bbecc78e",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, format='bshc', units='m')


__all__ = ['surface', 'bedrock', 'tbi', 'ret', 'ice']
