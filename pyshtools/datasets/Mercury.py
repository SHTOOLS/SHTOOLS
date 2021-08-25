'''
Datasets related to the planet Mercury.

Topography
----------
GTMES150           :  Topography constructed using laser altimeter and
                      occultation data.

Gravity
-------
JGMESS160A         :  Konolpiv et al. (2020)
JGMESS160A_ACCEL   :  Konolpiv et al. (2020)
JGMESS160A_TOPOSIG :  Konolpiv et al. (2020)
GGMES100V08        :  Genova et al. (2019)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ..shclasses import SHCoeffs as _SHCoeffs
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Mercury import omega as _omega


def GTMES150(lmax=150):
    '''
    GTMES150 is a GSFC 150 degree and order spherical harmonic model of the
    shape of the planet Mercury. This model is based on 26 million laser
    altimeter measurements and 557 radio occultations. The coefficients are in
    units of meters.

    Documentation for this model can be found on the PDS web site at
    https://pds-geosciences.wustl.edu/

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/gtmes_150v05_sha.tab",  # noqa: E501
        known_hash="sha256:c49d07c14b09c1b1ed1b4bfc4b42d2ff058875f5e949b50128b83ffa94c659b3",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    temp = _SHCoeffs.from_file(
        fname, lmax=lmax, header=True, errors=True, name='GTMES150', units='m',
        encoding='utf-8')
    temp.coeffs *= 1000.
    temp.errors *= 1000.
    return temp


def JGMESS160A(lmax=160):
    '''
    JGMESS160A is a JPL 160 degree and order spherical harmonic model of the
    gravitational potential of Mercury. This model applies a Kaula law
    constraint to all degrees.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Ermakov, A.I. (2020). The Mercury gravity
        field, orientation, Love number, and ephemeris from the MESSENGER
        radiometric tracking data, Icarus, 335, 253-260,
        doi:10.1016/j.icarus.2019.07.020.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.tab",  # noqa: E501
        known_hash="sha256:14fa0129c4b5ef655e08a883a05a476a836a806349da607f84b3c2b2e3d899ca",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='JGMESS160A', encoding='utf-8')


def JGMESS160A_ACCEL(lmax=160):
    '''
    JGMESS160A_ACCEL is a JPL 160 degree and order spherical harmonic model of
    the gravitational potential of Mercury. This model applies a surface
    acceleration constraint that is based on the degree strength given by the
    coefficient covariance.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Ermakov, A.I. (2020). The Mercury gravity
        field, orientation, Love number, and ephemeris from the MESSENGER
        radiometric tracking data, Icarus, 335, 253-260,
        doi:10.1016/j.icarus.2019.07.020.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_accel_sha.tab",  # noqa: E501
        known_hash="sha256:a0f99553fcea3d7d0c1395c89d8516f4e53e9e39893a2b82541ca1520b291423",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='JGMESS160A_ACCEL', encoding='utf-8')


def JGMESS160A_TOPOSIG(lmax=160):
    '''
    JGMESS160A_TOPOSIG is a JPL 160 degree and order spherical harmonic model
    of the gravitational potential of Mercury. This model applies a constraint
    similar to the Kaula constraint except that the uncertainty is given by the
    corresponding magnitude of the gravity derived from topography coefficient.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Ermakov, A.I. (2020). The Mercury gravity
        field, orientation, Love number, and ephemeris from the MESSENGER
        radiometric tracking data, Icarus, 335, 253-260,
        doi:10.1016/j.icarus.2019.07.020.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_toposig_sha.tab",  # noqa: E501
        known_hash="sha256:ad6d77e55968b9ddaea9cde03cb00ad6b630c81e8509e413e5b5b11213b3d848",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='JGMESS160A_TOPOSIG', encoding='utf-8')


def GGMES100V08(lmax=100):
    '''
    GGMES100V08 is a GSFC 100 degree and order spherical harmonic model of
    the gravitational potential of Mercury. This model applies a Kaula
    constraint to all degrees.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Genova, A., Goossens, S., Mazarico, E., Lemoine, F.G., Neumann, G.A.,
        Kuang, W., Sabaka, T.J., Hauck, S.A. II, Smith, D.E., Solomon S.C.,
        Zuber, M.T. (2019). Geodetic evidence that Mercury has a solid inner
        core. Geophysical Research Letters, 46, 3625-3633,
        doi:10.1029/2018GL081135.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/ggmes_100v08_sha.tab",  # noqa: E501
        known_hash="sha256:f81b33663ced0c6e05c775aa2d8c6c8eb99c41a20d4063d240fed536e45058fd",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='GGMES100V08', encoding='utf-8')


__all__ = ['GTMES150', 'JGMESS160A', 'JGMESS160A_ACCEL', 'JGMESS160A_TOPOSIG',
           'GGMES100V08']
