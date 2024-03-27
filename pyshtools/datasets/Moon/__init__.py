'''
Datasets related to Earth's Moon.

Shape
-----
LOLA_shape_pa     :  Wieczorek (2024)
LOLA_shape        :  Wieczorek (2024)

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
from pooch import create as _create
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import DOIDownloader as _DOIDownloader
from ...shclasses import SHCoeffs as _SHCoeffs
from ...shclasses import SHGravCoeffs as _SHGravCoeffs
from ...shclasses import SHMagCoeffs as _SHMagCoeffs
from ...constants.Moon import angular_velocity as _omega
from . import historical  # noqa: F401


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


def LOLA_shape(lmax=719):
    '''
    LOLA_shape is a spherical harmonic model of the shape of Earth's Moon in
    the mean Earth/polar axis coordinate system based on LOLA laser altimetry
    data obtained by the Lunar Reconaissance Orbiter mission. The maximum
    spherical harmonic degree of the model is 5759, which has an effective
    spatial resolution of 64 pixels per degree. Three lower resolution models
    are available in this archive (with lmax of 719, 1439 and 2879), and only
    the smallest that is required by the user input lmax will be downloaded. If
    lmax is not specified, the lowest resolution model (719) will be returned.
    If a negative value for lmax is specified, the maximum resolution model
    will be returned. The coefficients are in units of meters.

    This shape model should not be used in conjuction with most lunar gravity
    models, which use a principal axis coordinate system. For a principal axis
    model, use LOLA_shape_pa instead.

    Parameters
    ----------
    lmax : int, optional, default = 719
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Wieczorek, M. (2024). Spherical harmonic models of the shape of the Moon
        [LOLA] (1.0.1) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.10820774
    LRO LOLA Team (2013). LRO-L-LOLA-4-GDR-V1.0, NASA Planetary Data System.
    '''
    archive = _create(
        path=_os_cache('pyshtools'),
        base_url="doi:10.5281/zenodo.10820774",
        registry={
            "Moon_LOLA_shape_5759.bshc.gz": "sha256:7f1eddfccd007c56e983c43c8dc470fb58f8b23ede655f5bf4bf2208635c66ca",  # noqa: E501
            "Moon_LOLA_shape_2879.bshc.gz": "sha256:b22c3c71e2c14bc84e85304471c4b583542ee4e8f54ea2501bdc53d86017fe1f",  # noqa: E501
            "Moon_LOLA_shape_1439.bshc.gz": "sha256:f5e8f1b7acaab626db5828f9dc6f9ea1ce948bb5000e5bfcf1ff69b3141f1845",  # noqa: E501
            "Moon_LOLA_shape_719.bshc.gz": "sha256:20cafa143fda3e5ba5af6868c9e95054fabf68d0050c4be64a83d44df1333b4c",  # noqa: E501
            },
        )

    if lmax < 0:
        lmax = 5759

    if lmax >= 0 and lmax <= 719:
        fname = archive.fetch("Moon_LOLA_shape_719.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 719 and lmax <= 1439:
        fname = archive.fetch("Moon_LOLA_shape_1439.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 1439 and lmax <= 2879:
        fname = archive.fetch("Moon_LOLA_shape_2879.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
    else:
        fname = archive.fetch("Moon_LOLA_shape_5759.bshc.gz",
                              downloader=_DOIDownloader(progressbar=True))
        lmax = min(lmax, 5759)

    return _SHCoeffs.from_file(fname, lmax=lmax, name='LOLA_shape (Moon)',
                               units='m', format='bshc')


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
        url="doi:10.5281/zenodo.3873648/T2015_449.sh.gz",
        known_hash="sha256:4db0b77b3863f38d6fb6e62c5c1116bf7123b77c5aad65df7dae598714edd655",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHMagCoeffs.from_file(fname, lmax=lmax, header=True,
                                  file_units='T', name='T2015_449 (Moon)',
                                  units='nT', encoding='utf-8')


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
                                  header_units='km', name='Ravat2020 (Moon)',
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
                                   name='GRGM900C (Moon)', encoding='utf-8')


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
        url="https://pgda.gsfc.nasa.gov/data/MoonRM1/sha.grgm1200b_sigma",
        known_hash="sha256:f08a988b43f3eaa5a2089045a9b7e41e02f16542c7912b87ea34366fafa39bc5",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=1, gm_index=0, errors=True,
                                   omega=_omega.value, name='GRGM1200B (Moon)',
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
                                   name='GRGM1200B_RM1_1E0 (Moon)',
                                   encoding='utf-8')


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
                                   name='GL0900D (Moon)', encoding='utf-8')


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
                                   name='GL1500E (Moon)', encoding='utf-8')


__all__ = ['LOLA_shape_pa', 'LOLA_shape', 'T2015_449', 'Ravat2020', 'GRGM900C',
           'GRGM1200B', 'GRGM1200B_RM1_1E0', 'GL0900D', 'GL1500E',
           'historical']
