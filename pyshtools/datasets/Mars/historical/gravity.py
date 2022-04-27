'''
Historical datasets related to the gravity field of Mars.

Gravity
-------
GGMRO95A         :  Lemoine et al. (2008)
MRO95A           :  Konopliv et al. (2011)
MRO110B          :  Konopliv et al. (2011)
MRO110B2         :  Konopliv et al. (2011)
MRO110C          :  Konopliv et al. (2011)
MRO120D          :  Konopliv et al. (2016)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ....shclasses import SHGravCoeffs as _SHGravCoeffs
from ....constants.Mars import omega as _omega


def GGMRO95A(lmax=95):
    '''
    GGMRO95A is a GSFC 95 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint for
    degrees greater than 50.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Lemoine, F. G., Mazarico, E., Neumann, G., and D. Chinn (2008). New
        Solutions for the Mars Static and Temporal Gravity Field using the
        Mars Reconnaissance Orbiter Eos Trans. AGU, 89(53), Fall Meet. Suppl.,
        Abstract P41B-1376.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/ggmro_095a_sha.tab",  # noqa: E501
        known_hash="sha256:94978102bc6f443ff195dda0f9020660994381b1b5b47d9ad3110922e87b7ebc",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    clm = _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                  r0_index=0, gm_index=1, errors=True,
                                  omega=_omega.value, name='GGMRO95A',
                                  encoding='utf-8')
    clm.r0 *= 1.e3
    return clm


def MRO95A(lmax=95):
    '''
    MRO95A is a JPL 95 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint for
    degrees greater than 70.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Asmar, S.W., Folkner, W.M., Karatekin, O., Nunes, D.C.,
        Smrekar, S.E., Yoder, C.F., Zuber, M.T. (2011). Mars high resolution
        gravity fields from MRO, Mars seasonal gravity, and other dynamical
        parameters. Icarus, 211, 401-428, doi:10.1016/j.icarus.2010.10.004.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_095a_sha.tab",  # noqa: E501
        known_hash="sha256:c35927f287f40bd039d101e50b5364841ebfd9ab55d0a863784175c5bc984ce5",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   omega=_omega.value, name='MRO95A',
                                   encoding='utf-8')


def MRO110B(lmax=110):
    '''
    MRO110B is a JPL 110 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint for
    degrees greater than 70.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Asmar, S.W., Folkner, W.M., Karatekin, O., Nunes, D.C.,
        Smrekar, S.E., Yoder, C.F., Zuber, M.T. (2011). Mars high resolution
        gravity fields from MRO, Mars seasonal gravity, and other dynamical
        parameters. Icarus, 211, 401-428, doi:10.1016/j.icarus.2010.10.004.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_110b_sha.tab",  # noqa: E501
        known_hash="sha256:49e55b6833ed90c7853abbdb0b54bcc30a3613b20cc7dbbfcb7399a0e3b70e5c",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   omega=_omega.value, name='MRO110B',
                                   encoding='utf-8')


def MRO110B2(lmax=110):
    '''
    MRO110B2 is a JPL 110 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint that
    is 35% looser than used in MRO110B for degrees greater than 80 (as opposed
    to 70 for MRO110B). The data included in this solution are same as used in
    model MRO110B.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Asmar, S.W., Folkner, W.M., Karatekin, O., Nunes, D.C.,
        Smrekar, S.E., Yoder, C.F., Zuber, M.T. (2011). Mars high resolution
        gravity fields from MRO, Mars seasonal gravity, and other dynamical
        parameters. Icarus, 211, 401-428, doi:10.1016/j.icarus.2010.10.004.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_110b2_sha.tab",  # noqa: E501
        known_hash="sha256:e1fc6dc9578715a328a43407d8c52f173a1c72740981a46e107530e05c813d1c",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   omega=_omega.value, name='MRO110B2',
                                   encoding='utf-8')


def MRO110C(lmax=110):
    '''
    MRO110C is a JPL 110 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint for
    degrees greater than 70.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Asmar, S.W., Folkner, W.M., Karatekin, O., Nunes, D.C.,
        Smrekar, S.E., Yoder, C.F., Zuber, M.T. (2011). Mars high resolution
        gravity fields from MRO, Mars seasonal gravity, and other dynamical
        parameters. Icarus, 211, 401-428, doi:10.1016/j.icarus.2010.10.004.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_110c_sha.tab",  # noqa: E501
        known_hash="sha256:31fecb73a8959a55a9e5f92911d23ddaecb0a04139b93b500406b7944a6d6b80",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   omega=_omega.value, name='MRO110C',
                                   encoding='utf-8')


def MRO120D(lmax=120):
    '''
    MRO120D is a JPL 120 degree and order spherical harmonic model of the
    gravitational potential of Mars. This model applies a Kaula constraint for
    degrees greater than 80.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Folkner, W.M. (2016). An improved JPL Mars
        gravity field and orientation from Mars orbiter and lander tracking
        data. Icarus, 274, 253-260, doi:10.1016/j.icarus.2016.02.052.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120d_sha.tab",  # noqa: E501
        known_hash="sha256:00c3a2fada7bdfb8022962752b1226be1de21ea469f9e88af5d8dc47d23883bd",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   omega=_omega.value, name='MRO120D',
                                   encoding='utf-8')


__all__ = ['GGMRO95A', 'MRO95A', 'MRO110B', 'MRO110B2', 'MRO110C', 'MRO120D']
