'''
Historical datasets related to the shape of Venus.

Topography
-------
SHTJV360A01 (Magellan)       :  Rappaport and Plaut (1994)
SHTJV360A02 (Magellan)       :  Rappaport et al. (1999)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ....shclasses import SHCoeffs as _SHCoeffs


def SHTJV360A01(lmax=360):
    '''
    SHTJV360A01 is a 360 degree and order spherical harmonic model of the
    shape of the planet Venus derived exclusively from Magellan altimetry.
    The coefficients are in units of meters.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Rappaport, N. J., Plaut, J. J. (1994). A 360 degree and order model of
        Venus topography, Icarus, 112, 27-33, doi:10.1006/icar.1994.1167.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/topo/shtjv360.a01",  # noqa: E501
        known_hash="sha256:eafc429a9a13b11d4142f19af2092d40dc7a69fb4ee915d99b5738f9b08171cd",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    clm = _SHCoeffs.from_file(fname, lmax=lmax, name='SHTJV360A01',
                              units='m', encoding='utf-8', skip=240,
                              header=True)
    clm *= float(clm.header[1]) * 1.e3
    return clm


def SHTJV360A02(lmax=360):
    '''
    SHTJV360A02 is a 360 degree and order spherical harmonic model of the
    shape of the planet Venus derived exclusively from Magellan altimetry.
    The coefficients are in units of meters. This version is an update of the
    previous model SHTJV360A01 that used more accurate Magellan spacecraft
    ephemerides.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Rappaport, N. J., Konopliv, A. S., and Kucinskas, A. B. (1999). An improved
        360 degree and order model of Venus topography, Icarus, 139, 19-31,
        doi:10.1006/icar.1999.6081.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/topo/shtjv360.a02",  # noqa: E501
        known_hash="sha256:c5cc0da9890826b0531f57ed408a64be9fe3a679c81a3e8f6405ad5c8fdcdeb4",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    clm = _SHCoeffs.from_file(fname, lmax=lmax, name='SHTJV360A02',
                              units='m', encoding='utf-8', skip=242,
                              header=True)
    rt = float(clm.header[1]) * 1.e3
    clm *= rt
    clm.coeffs[0, 0, 0] = rt
    return clm


__all__ = ['SHTJV360A01', 'SHTJV360A02']
