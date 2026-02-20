'''
Datasets related to the asteroid (4) Vesta.

Shape
-----
DLR_SPG_shape  :  Wieczorek (2024)

Gravity
-------
VESTA20H       :  Konopliv et al. (2014)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import create as _create
from pooch import HTTPDownloader as _HTTPDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..shclasses import SHCoeffs as _SHCoeffs
from ..constants.Vesta import angular_velocity as _omega
from ._utils import _choose_sh_model


def DLR_SPG_shape(lmax=719):
    '''
    DLR_SPG_shape is a spherical harmonic model of the shape of asteroid (4)
    Vesta based on stereo photogrammetric data obtained by the Dawn mission.
    The maximum spherical harmonic degree of the model is 5759, which has an
    effective spatial resolution of 64 pixels per degree. Three lower
    resolution models are available in this archive (with lmax of 719, 1439
    and 2879), and only the smallest that is required by the user input lmax
    will be downloaded. If lmax is not specified, the lowest resolution model
    (719) will be returned. If a negative value for lmax is specified, the
    maximum resolution model will be returned. The coefficients are in units
    of meters. This model uses the IAU Claudia double-prime coordinate system.

    Parameters
    ----------
    lmax : int, optional, default = 719
        The maximum spherical harmonic degree to return.

    References
    ----------
    Wieczorek, M. (2024). Spherical harmonic models of the shape of the
        asteroid (4) Vesta [DLR SPG] (1.0.2) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.10820681
    Preusker, F., F. Scholten, K.-D Matz, T. Roatsch, R. Jaumann, C.A. Raymond,
        and C.T. Russell (2016). DAWN FC2 DERIVED VESTA DTM SPG V1.0,
        DAWN-A-FC2-5-VESTADTMSPG-V1.0, NASA Planetary Data System.
    '''
    archive = _create(
        path=_os_cache('pyshtools'),
        base_url="doi:10.5281/zenodo.10820681",
        registry={
            "Vesta_DLR_SPG_shape_5759.bshc.gz": "sha256:3a2486575fc99469bfe243b78c92cc3c0b968cbf1f278436079ef53961a75670",  # noqa: E501
            "Vesta_DLR_SPG_shape_2879.bshc.gz": "sha256:04deabadc25338cc4de5320447afbcc429baa5715427528641d8e0561521f657",  # noqa: E501
            "Vesta_DLR_SPG_shape_1439.bshc.gz": "sha256:99eae34532e93a4611be51436de7e3e9be8c3f64a799bc43444013ecf73ca6d3",  # noqa: E501
            "Vesta_DLR_SPG_shape_719.bshc.gz": "sha256:140532d5ca7070e677ac064b9ff59c60baa21811b0541138e564e94ec8eafec9",  # noqa: E501
            },
        )

    fname, lmax = _choose_sh_model(
        archive=archive,
        user_lmax=lmax,
    )

    return _SHCoeffs.from_file(fname, lmax=lmax, name='DLR_SPG_shape (Vesta)',
                               units='m', format='bshc')


def VESTA20H(lmax=20):
    '''
    VESTA20H is a JPL 20 degree and order spherical harmonic model of the
    gravitational potential of asteroid (4) Vesta. This model uses a coordinate
    system that is similar to the IAU Claudia double-prime coordinate system.

    Parameters
    ----------
    lmax : int, optional, default = 20
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Asmar, S.W., Park, R.S., Bills, B.G., Centinello, F.,
        Chamberlin, A.B., Ermakov, A., Gaskell, R.W., Rambaux, N.,
        Raymond, C.A., Russell, C.T., Smith, D.E., Tricarico, P., Zuber, M.T.
        (2014). The Vesta gravity field, spin pole and rotation period,
        landmark positions, and ephemeris from the Dawn tracking and optical
        data. Icarus, 240, 103-117, doi:10.1016/j.icarus.2013.09.005.
    '''
    if lmax < 0:
        lmax = 20

    fname = _retrieve(
        url="https://sbnarchive.psi.edu/pds3/dawn/grav/DWNVGRS_2/DATA/SHADR/JGDWN_VES20H_SHA.TAB",  # noqa: E501
        known_hash="sha256:389bb30882b43b36877e9ff8bf54c55a210d7ba98cab88e9c4a1fa83ad8484b7",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='VESTA20H (Vesta)', encoding='utf-8',
                                   omega=_omega.value)


__all__ = ['DLR_SPG_shape', 'VESTA20H']
