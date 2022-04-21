'''
Datasets related to the planet Mars.

Topography
----------
MarsTopo2600     :  Wieczorek (2015)
MarsTopo719      :  Wieczorek (2015)

Gravity
-------
GMM3             :  Genova et al. (2016)
GMM3_RM1_1E0     :  Goossens et al. (2017)
MRO120F          :  Konopliv et al. (2020)

Magnetic field
--------------
Langlais2019     :  Langlais et al. (2019)
Morschhauser2014 :  Morschhauser et al. (2014)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ...shclasses import SHCoeffs as _SHCoeffs
from ...shclasses import SHGravCoeffs as _SHGravCoeffs
from ...shclasses import SHMagCoeffs as _SHMagCoeffs
from pooch import Decompress as _Decompress
from ...constants.Mars import omega as _omega
from . import historical  # noqa: F401


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
        url="https://zenodo.org/record/3870922/files/MarsTopo2600.shape.gz",
        known_hash="sha256:8882a9ee7ee405d971b752028409f69bd934ba5411f1c64eaacd149e3b8642af",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
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
        url="https://zenodo.org/record/6475460/files/MarsTopo719.shape.gz",
        known_hash="sha256:37a98efae5eab7c85260f4b43315fe9fcf44247a61581bed1b6f7f10f79adea0",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='MarsTopo719',
                               units='m', encoding='utf-8')


def GMM3(lmax=120):
    '''
    GMM3 is a GSFC 120 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint for
    degrees greater than 90.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Genova, A., Goossens, S., Lemoine, F.G., Mazarico, E., Neumann, G.A.,
        Smith, D.E., Zuber, M.T. (2016). Seasonal and static gravity field of
        Mars from MGS, Mars Odyssey and MRO radio science, Icarus, 272,
        228-245, doi:10.1016/j.icarus.2016.02.050.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/gmm3_120_sha.tab",  # noqa: E501
        known_hash="sha256:eb4913b1afb6682406e6a9dad5be7918a162fa8462473c9a2e7aae258d4c2c9c",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   omega=_omega.value, name='GMM3',
                                   encoding='utf-8')


def GMM3_RM1_1E0(lmax=150):
    '''
    GMM3_RM1_1E0 is a GSFC 150 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model uses the same data as GMM3, but
    with a rank-minus-1 constraint based on gravity from surface topography for
    degrees greater than 50 with a value of lambda equal to 1.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Goossens, S., Sabaka, T.J., Genova, A, Mazarico, E., Nicholas, J.B.,
        Neumann, G.A. (2017), Evidence for a low bulk crustal density for Mars
        from gravity and topography, Geophysical Research Letters, 44,
        7686-7694, doi:10.1002/2017GL074172.
    '''
    fname = _retrieve(
        url="https://pgda.gsfc.nasa.gov/data/MarsDensityRM1/sha.gmm3_l150_rm1_lambda_1",  # noqa: E501
        known_hash="sha256:b309917362bd2014df42a62cb19ea321ee8db97997b0688eda2774deb46ef538",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=1, gm_index=0, errors=False,
                                   omega=_omega.value, name='GMM3_RM1_1E0',
                                   encoding='utf-8')


def MRO120F(lmax=120):
    '''
    MRO120F is a JPL 120 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint for
    degrees greater than 80.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Rivoldini, A., Baland, R.-M., Le Maistre, S.,
        Van Hoolst, T., Yseboodt, M., Dehant, V. (2020). Detection of the
        Chandler Wobble of Mars From Orbiting Spacecraft. Geophysical Research
        Letters, 47, e2020GL090568, doi:10.1029/2020GL090568.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120f_sha.tab",  # noqa: E501
        known_hash="sha256:ddd3de9c30d75879fe37aa17a1149e7c96c141c095962954cb7ea865a2c025b6",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   omega=_omega.value, name='MRO120F',
                                   encoding='utf-8')


def Langlais2019(lmax=134):
    '''
    Langlais2019 is a 134 degree and order spherical harmonic model of the
    magnetic potential of Mars. This model makes use of data from MGS MAG,
    MGS ER and MAVEN MAG. The coefficients are output in units of nT.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    References
    ----------
    Langlais, B., Thébault, E., Houliez, A., Purucker, M.E., Lillis, R.J.
        (2019). A new model of the crustal magnetic field of Mars using MGS
        and MAVEN. Journal of Geophysical Research: Planets, 124, 1542-1569,
        doi:10.1029/2018JE005854.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/3876714/files/Langlais2019.sh.gz?download=1",  # noqa: E501
        known_hash="sha256:3cad9e268f0673be1702f1df504a4cbcb8dba4480c7b3f629921911488fe247b",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHMagCoeffs.from_file(fname, lmax=lmax, skip=4, r0=3393.5e3,
                                  header=False, file_units='nT',
                                  name='Langlais2019', units='nT',
                                  encoding='utf-8')


def Morschhauser2014(lmax=110):
    '''
    Morschhauser2014 is a 110 degree and order spherical harmonic model of the
    magnetic potential of Mars. The coefficients are output in units of nT.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    References
    ----------
    Morschhauser, A., Lesur, V., Grott, M. (2014). A spherical harmonic model
        of the lithospheric magnetic ﬁeld of Mars, Journal of Geophysical
        Research: Planets, 119, 1162-1188, doi:10.1002/2013JE004555.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/3876495/files/Morschhauser2014.txt.gz?download=1",  # noqa: E501
        known_hash="sha256:a86200b3147a24447ff8bba88ec6047329823275813a9f5e9505bb611e3e86e0",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
        processor=_Decompress(),
    )
    return _SHMagCoeffs.from_file(fname, r0=3393.5e3, skip=3, header=False,
                                  format='dov', file_units='nT',
                                  name='Morschhauser2014', units='nT',
                                  encoding='utf-8')


__all__ = ['MarsTopo2600', 'GMM3', 'GMM3_RM1_1E0', 'MRO120F', 'Langlais2019',
           'Morschhauser2014', 'historical', 'MarsTopo719']
