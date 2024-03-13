'''
Datasets related to the dwarf planet (1) Ceres.

Shape
-----
DLR_SPG_shape  :  Wieczorek (2024)
JPL_SPC_shape  :  Wieczorek (2024)

Gravity
-------
CERES18D       :  Konopliv et al. (2018)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import create as _create
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..shclasses import SHCoeffs as _SHCoeffs
from ..constants.Ceres import omega as _omega


def DLR_SPG_shape(lmax=719):
    '''
    DLR_SPG_shape is a spherical harmonic model of the shape of asteroid (1)
    Ceres based on stereo photogrammetric data obtained by the Dawn mission.
    The maximum spherical harmonic degree of the model is 5399, which has an
    effective spatial resolution of 60 pixels per degree. Three lower
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

    References
    ----------
    Wieczorek, M. (2024). Spherical harmonic models of the shape of the
        asteroid (1) Ceres [DLR SPG] (1.0.0) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.10804157
    Roatsch, T., E. Kersten, K.-D. Matz, F. Preusker, F. Scholten, S. Elgner,
        S.E. Schroeder, R. Jaumann, C.A. Raymond, C.T. Russell (2016). DAWN FC2
        DERIVED CERES HAMO DTM SPG V1.0, DAWN-A-FC2-5-CERESHAMODTMSPG-V1.0,
        NASA Planetary Data System
    '''
    archive = _create(
        path=_os_cache('pyshtools'),
        base_url="doi:10.5281/zenodo.10804157",
        registry={
            "Ceres_shape_5399.sh.gz": "sha256:a474a8ccc7620673c1516444a815762bbd3cd37f7be4fbe0634cb781c18c355f",  # noqa: E501
            "Ceres_shape_2879.sh.gz": "sha256:9ec1442275fd6b0a19830766cfe023c672133b060850324cfc36fd445dcf0f9e",  # noqa: E501
            "Ceres_shape_1439.sh.gz": "sha256:ae3e99187ccfae9c986772f4d0f9e932d51ea292f9d60214589468aae12d8f8e",  # noqa: E501
            "Ceres_shape_719.sh.gz": "sha256:e3dbbfb7c1cc55b788f2a8ca751f0df727880f463fabfd7ec6ac9cb15b9a2b96",  # noqa: E501
            },
        )

    if lmax < 0:
        lmax = 5399

    if lmax >= 0 and lmax <= 719:
        fname = archive.fetch("Ceres_shape_719.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 719 and lmax <= 1439:
        fname = archive.fetch("Ceres_shape_1439.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    elif lmax > 1439 and lmax <= 2879:
        fname = archive.fetch("Ceres_shape_2879.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    else:
        fname = archive.fetch("Ceres_shape_5399.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
        lmax = min(lmax, 5399)

    return _SHCoeffs.from_file(fname, lmax=lmax, name='DLR_SPG_shape (Ceres)',
                               units='m', format='bshc')


def JPL_SPC_shape(lmax=719):
    '''
    JPL_SPC_shape is a spherical harmonic model of the shape of asteroid (1)
    Ceres based on stereo photoclinometric data obtained by the Dawn mission.
    The maximum spherical harmonic degree of the model is 1023, which has an
    effective spatial resolution of 11.3 pixels per degree. One lower
    resolution model is available in this archive with lmax of 719, and only
    the smallest that is required by the user input lmax will be downloaded.
    If lmax is not specified, the lowest resolution model (719) will be
    returned. If a negative value for lmax is specified, the maximum resolution
    model will be returned. The coefficients are in units of meters.

    Parameters
    ----------
    lmax : int, optional, default = 719
        The maximum spherical harmonic degree to return.

    References
    ----------
    JPL_SPC_shape: Wieczorek, M. (2024). Spherical harmonic models of the
        shape of asteroid (1) Ceres [JPL SPC] (1.0.0) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.10812848
    Park, R.S. and Buccino, D.R., Ceres SPC Shape Model Dataset V1.0.
        DAWN-A-FC2-5-CERESSHAPESPC-V1.0. NASA Planetary Data System, 2018
    Park, R. S., Vaughan, A. T., Konopliv, A. S., Ermakov, A. I., Mastrodemos,
        N., Castillo-Rogez, J. C., Joy, S. P., Nathues, A., Polanskey, C. A.,
        Rayman, M. D., Riedel, J. E., Raymond, C. A., Russell, C. T., & Zuber,
        M. T. (2019). High-resolution shape model of Ceres from
        stereophotoclinometry using Dawn Imaging Data. Icarus, 319, 812–827.
        https://doi.org/10.1016/j.icarus.2018.10.024
    '''
    archive = _create(
        path=_os_cache('pyshtools'),
        base_url="doi:10.5281/zenodo.10812848",
        registry={
            "Ceres_JPL_SPC_shape_1023.sh.gz": "sha256:9ca1b3c31760beba01c56ac7f2c1d30d62b1480aad551eebf082f0e34eb19f06",  # noqa: E501
            "Ceres_JPL_SPC_shape_719.sh.gz": "sha256:5e66eeeb96bfbdfc30e8de7f13f6a48dd8795c7a3f781786bb985bf5c47572b5",  # noqa: E501
            },
        )

    if lmax < 0:
        lmax = 1023

    if lmax >= 0 and lmax <= 719:
        fname = archive.fetch("Ceres_JPL_SPC_shape_719.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    else:
        fname = archive.fetch("Ceres_JPL_SPC_shape_1023.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
        lmax = min(lmax, 1023)

    return _SHCoeffs.from_file(fname, lmax=lmax, name='JPL_SPC_shape (Ceres)',
                               units='m', format='bshc')


def CERES18D(lmax=18):
    '''
    CERES18D is a JPL 18 degree and order spherical harmonic model of the
    gravitational potential of (1) Ceres.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Vaughan, A.T., Bills, B.G., Asmar, S.W.,
        Ermakov, A.I., Rambaux, N., Raymond, C.A., Castillo-Rogez, J.C.,
        Russell, C.T., Smith, D.E., Zuber, M.T. (2018). The Ceres gravity
        field, spin pole, rotation period and orbit from the Dawn radiometric
        tracking and optical data, Icarus, 299, 411-429,
        doi:10.1016/j.icarus.2017.08.005.
    '''
    fname = _retrieve(
        url="https://sbnarchive.psi.edu/pds3/dawn/grav/DWNCGRS_2_v3_181005/DATA/SHADR/JGDWN_CER18D_SHA.TAB",  # noqa: E501
        known_hash="sha256:e7ccb1f0c689f77fe5dae4e0bb5d514db1cf5acb5be927bbcaa8576aca153981",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='CERES18D', encoding='utf-8',
                                   omega=_omega.value)


__all__ = ['DLR_SPG_shape', 'CERES18D']
