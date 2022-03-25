'''
Datasets related to Earth's Moon.

Topography
----------
MoonTopo2600p     :  Wieczorek (2015)

Gravity
-------
GRGM900C          :  Lemoine et al. (2014)
GRGM1200B         :  Goossens et al. (2020)
GRGM1200B_RM1_1E0 :  Goossens et al. (2020)
GL0900D           :  Konopliv et al. (2014)
GL1500E           :  Konopliv et al. (2014)

Magnetic field
--------------
T2015_449         :  Wieczorek (2018), using data from Tsunakawa et al. (2015)
Ravat2020         :  Ravat et al. (2020)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ...shclasses import SHCoeffs as _SHCoeffs
from ...shclasses import SHGravCoeffs as _SHGravCoeffs
from ...shclasses import SHMagCoeffs as _SHMagCoeffs
from ...constants.Moon import omega as _omega
from . import historical  # noqa: F401


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
        url="https://zenodo.org/record/3870924/files/MoonTopo2600p.shape.gz",
        known_hash="sha256:193146df894e2fef796df9d6142c78fae6fa5c183fd79d3f79eeb356602af69a",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='MoonTopo2600p',
                               units='m', encoding='utf-8')


def T2015_449(lmax=449):
    '''
    T2015_449 is a 449 degree and order spherical harmonic model of the
    magnetic potential of the Moon. This model was used in Wieczorek (2018) and
    is a spherical harmonic expansion of the global magnetic field model of
    Tsunakawa et al. (2015). The coefficients are output in units of nT.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    References
    ----------
    Tsunakawa, H., Takahashi, F., Shimizu, H., Shibuya, H., Matsushima, M.
        (2015). Surface vector mapping of magnetic anomalies over the Moon
        using Kaguya and Lunar Prospector observations, Journal of Geophysical
        Research Planets, 120, 1160-1185, doi:10.1002/2014JE004785.
    Wieczorek, M.A. (2018). Strength, depth, and geometry of magnetic sources
        in the crust of the Moon from localized power spectrum analysis,
        Journal of Geophysical Research Planets, 123, 291-316,
        doi:10.1002/2017JE005418.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/3873648/files/T2015_449.sh.gz",
        known_hash="sha256:4db0b77b3863f38d6fb6e62c5c1116bf7123b77c5aad65df7dae598714edd655",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHMagCoeffs.from_file(fname, lmax=lmax, header=True,
                                  file_units='T', name='T2015_449', units='nT',
                                  encoding='utf-8')


def Ravat2020(lmax=450):
    '''
    Ravat2020 is a 450 degree and order spherical harmonic model of the
    magnetic potential of the Moon. This model is based on using magnetic
    monopoles with an L1-norm regularisation. The coefficients are output
    in units of nT.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    References
    ----------
    Ravat, D., Purucker, M. E., Olsen, N. (2020). Lunar magnetic field models
        from Lunar Prospector and SELENE/Kaguya along‐track magnetic field
        gradients, Journal of Geophysical Research: Planets, 125,
        e2019JE006187, doi:10.1029/2019JE006187.
    '''
    fname = _retrieve(
        url="https://uknowledge.uky.edu/cgi/viewcontent.cgi?filename=4&article=1001&context=ees_data&type=additional",  # noqa: E501
        known_hash="sha256:dd1128d7819a8de097f3abeba93fee4cb80fced5bd63d56cca5a9bc70ac2bea9",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHMagCoeffs.from_file(fname, lmax=lmax, header=True, skip=8,
                                  header_units='km', name='Ravat2020',
                                  encoding='utf-8')


def GRGM900C(lmax=900):
    '''
    GRGM900C is a GSFC 900 degree and order spherical harmonic model of the
    gravitational potential of the Moon. This model applies a Kaula constraint
    for degrees greater than 600.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Lemoine, F.G., Goossens, S., Sabaka, T.J., Nicholas, J.B., Mazarico, E.,
        Rowlands, D.D., Loomis, B.D., Chinn, D.S., Neumann, G.A., Smith, D.E.,
        Zuber, M.T. (2014) GRGM900C: A degree 900 lunar gravity model from
        GRAIL primary and extended mission data, Geophysical Research Letters,
        41, 3382-3389, doi:10.1002/2014GL060027.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_0900c_sha.tab",  # noqa: E501
        known_hash="sha256:dab6ab06e0d3d7cbc594ea4bd03151a65534ed5fdf4f147ae38662428c04454e",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='GRGM900C', encoding='utf-8')


def GRGM1200B(lmax=1200):
    '''
    GRGM1200B is a GSFC 1200 degree and order spherical harmonic model of the
    gravitational potential of the Moon. This model applies a Kaula constraint
    for degrees greater than 600.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Goossens, S., Sabaka, T.J., Wieczorek, M.A., Neumann, G.A., Mazarico, E.,
        Lemoine, F.G., Nicholas, J.B., Smith. D.E., Zuber M.T. (2020).
        High-resolution gravity field models from GRAIL data and implications
        for models of the density structure of the Moon's crust. Journal of
        Geophysical Research Planets, 125, e2019JE006086,
        doi:10.1029/2019JE006086.
    '''
    fname = _retrieve(
        url="https://pgda.gsfc.nasa.gov/data/MoonRM1/sha.grgm1200b_sigma",  # noqa: E501
        known_hash="sha256:f08a988b43f3eaa5a2089045a9b7e41e02f16542c7912b87ea34366fafa39bc5",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=1, gm_index=0, errors=True,
                                   omega=_omega.value, name='GRGM1200B',
                                   encoding='utf-8')


def GRGM1200B_RM1_1E0(lmax=1200):
    '''
    GRGM1200B_RM1_1E0 is a GSFC 1200 degree and order spherical harmonic model
    of the gravitational potential of the Moon. This model uses a rank-minus-1
    constraint based on gravity from surface topography for degrees greater
    than 600 with a value of lambda equal to 1.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Goossens, S., Sabaka, T.J., Wieczorek, M.A., Neumann, G.A., Mazarico, E.,
        Lemoine, F.G., Nicholas, J.B., Smith. D.E., Zuber M.T. (2020).
        High-resolution gravity field models from GRAIL data and implications
        for models of the density structure of the Moon's crust. Journal of
        Geophysical Research Planets, 125, e2019JE006086,
        doi:10.1029/2019JE006086.
    '''
    fname = _retrieve(
        url="https://pgda.gsfc.nasa.gov/data/MoonRM1/sha.grgm1200b_rm1_1e0_sigma",  # noqa: E501
        known_hash="sha256:d42536cc716f5da8e067aa79a253c310e9d53d1d3b3ae7b43fa4517654d20d35",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=1, gm_index=0, errors=True,
                                   omega=_omega.value,
                                   name='GRGM1200B_RM1_1E0', encoding='utf-8')


def GL0900D(lmax=900):
    '''
    GL0900D is a JPL 900 degree and order spherical harmonic model of
    the gravitational potential of the Moon. This model applies a Kaula
    constraint for degrees greater than 700.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Yuan, D.-N., Asmar, S.W., Watkins, M.M.,
        Williams, J.G., Fahnestock, E., Kruizinga, G., Paik, M., Strekalov, D.,
        Harvey, N., Smith, D.E., Zuber, M.T. (2014). High-resolution lunar
        gravity fields from the GRAIL Primary and Extended Missions,
        Geophysical Research Letters, 41, 1452-1458, doi:10.1002/2013GL059066.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/jggrx_0900d_sha.tab",  # noqa: E501
        known_hash="sha256:0ead4e6260729c53fe29dcc6d954d473e87eb3ae1f4d932496a1417fd62fcdc6",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='GL0900D', encoding='utf-8')


def GL1500E(lmax=1500):
    '''
    GL1500E is a JPL 1500 degree and order spherical harmonic model of
    the gravitational potential of the Moon. This model applies a Kaula
    constraint for degrees greater than 700.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Yuan, D.-N., Asmar, S.W., Watkins, M.M.,
        Williams, J.G., Fahnestock, E., Kruizinga, G., Paik, M., Strekalov, D.,
        Harvey, N., Smith, D.E., Zuber, M.T. (2014). High-resolution lunar
        gravity fields from the GRAIL Primary and Extended Missions,
        Geophysical Research Letters, 41, 1452-1458, doi:10.1002/2013GL059066.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/jggrx_1500e_sha.tab",  # noqa: E501
        known_hash="sha256:93a7467b9241f6f94c131126c87fd3b81cfc3d223c474ee3444bf162b7c97f5a",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='GL1500E', encoding='utf-8')


__all__ = ['MoonTopo2600p', 'T2015_449', 'Ravat2020', 'GRGM900C', 'GRGM1200B',
           'GRGM1200B_RM1_1E0', 'GL0900D', 'GL1500E', 'historical']
