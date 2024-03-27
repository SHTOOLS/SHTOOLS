'''
Datasets related to the planet Mercury.

Shape
-----
USGS_SPG_shape     :  Maia (2024)
GTMES150           :  Shape model constructed using laser altimeter and
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
from pooch import create as _create
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHCoeffs as _SHCoeffs
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Mercury import angular_velocity as _omega


def USGS_SPG_shape(lmax=719):
    '''
    USGS_SPG_shape is a spherical harmonic model of the shape of Mercury based
    on stereo photogrammetric data obtained by the MESSENGER mission. The
    maximum spherical harmonic degree of the model is 5759, which has an
    effective spatial resolution of 64 pixels per degree. Three lower
    resolution models are available in this archive (with lmax of 719, 1439
    and 2879), and only the smallest that is required by the user input lmax
    will be downloaded. If lmax is not specified, the lowest resolution model
    (719) will be returned. If a negative value for lmax is specified, the
    maximum resolution model will be returned. The coefficients are in units
    of meters.

    Parameters
    ----------
    lmax : int, optional, default = 719
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Maia, J. (2024). Spherical harmonic models of the shape of Mercury
        [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10809345

    '''
    archive = _create(
        path=_os_cache('pyshtools'),
        base_url="doi:10.5281/zenodo.10809345",
        registry={
            "Mercury_shape_5759.sh.gz": "sha256:333ef33584cec5954933d1993c63c60cf2592298ed2300cb03fbac2037264a0b",  # noqa: E501
            "Mercury_shape_2879.sh.gz": "sha256:1e9dcfae93e8ba131bd73d1dc9fd31f36aca92004357be4785429cdb670ab5c1",  # noqa: E501
            "Mercury_shape_1439.sh.gz": "sha256:48d5ee575b5ffe97fd84b1c1dd2a0c1be0b59d476f52df23373bfd9810d59288",  # noqa: E501
            "Mercury_shape_719.sh.gz": "sha256:c481d89d3205dcda837868b25411f286d7f935a1a1f119a84558d88cb8b8da6b",  # noqa: E501
            },
        )

    if lmax < 0:
        lmax = 5759

    if lmax >= 0 and lmax <= 719:
        fname = archive.fetch("Mercury_shape_719.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 719 and lmax <= 1439:
        fname = archive.fetch("Mercury_shape_1439.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 1439 and lmax <= 2879:
        fname = archive.fetch("Mercury_shape_2879.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    else:
        fname = archive.fetch("Mercury_shape_5759.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
        lmax = min(lmax, 5759)

    return _SHCoeffs.from_file(fname, lmax=lmax,
                               name='USGS SPG_shape (Mercury)', units='m',
                               format='bshc')


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
        fname, lmax=lmax, header=True, errors=True, name='GTMES150 (Mercury)',
        units='m', encoding='utf-8')
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
                                   name='JGMESS160A (Mercury)',
                                   encoding='utf-8')


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
                                   name='JGMESS160A_ACCEL (Mercury)',
                                   encoding='utf-8')


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
                                   name='JGMESS160A_TOPOSIG (Mercury)',
                                   encoding='utf-8')


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
                                   name='GGMES100V08 (Mercury)',
                                   encoding='utf-8')


__all__ = ['USGS_SPG_shape', 'GTMES150', 'JGMESS160A', 'JGMESS160A_ACCEL',
           'JGMESS160A_TOPOSIG', 'GGMES100V08']
