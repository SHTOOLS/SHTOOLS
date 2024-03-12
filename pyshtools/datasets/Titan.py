'''
Datasets related to Saturn's moon Titan.

Shape
-----
Corlies2017_shape  :  Corlies et al. (2017)
Mitri2014_shape    :  Mitri et al. (2014)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHCoeffs as _SHCoeffs
from ..constants.Titan import omega as _omega


def Corlies2017_shape(lmax=8):
    '''
    Corlies2017_shape is a spherical harmonic model of the shape of Saturn's
    moon Titan based on radar altimeter, SAR topography, and radar stereo-
    photogrammetric data acquired by the Cassini mission. The maximum spherical
    harmonic degree of the model is 8, and the coefficients are in units of
    meters.

    Parameters
    ----------
    lmax : int, optional, default = 8
        The maximum spherical harmonic degree to return.

    References
    ---------
    Corlies, P., Hayes, A. G., Birch, S. P. D., Lorenz, R., Stiles, B. W.,
        Kirk, R., Poggiali, V., Zebker, H., & Iess, L. (2017). Titan’s
        Topography and Shape at the End of the Cassini Mission. Geophysical
        Research Letters, 44(23), 11,754-11,761.
        https://doi.org/10.1002/2017GL075518

    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10806604/Titan_shape_Corlies2017_unnorm.sh",  # noqa: E501
        known_hash="sha256:413c147243ea8e171c14161a0972dcfcf0980ce57039ab7d70714f7e2d1e94c0",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )

    return _SHCoeffs.from_file(fname, lmax=lmax, name='Corlies2017_shape',
                               units='m', errors=True, format='dov',
                               encoding='utf-8', normalization='unnorm',
                               csphase=-1)


def Mitri2014_shape(lmax=6):
    '''
    Mitri2014_shape is a spherical harmonic model of the shape of Saturn's moon
    Titan based on radar altimeter and SAR topography data acquired by the
    Cassini mission. The maximum spherical harmonic degree of the model is 6,
    and the coefficients are in units of meters.

    Parameters
    ----------
    lmax : int, optional, default = 6
        The maximum spherical harmonic degree to return.

    References
    ---------
    Mitri, G., Meriggiola, R., Hayes, A., Lefevre, A., Tobie, G., Genova, A.,
        Lunine, J. I., & Zebker, H. (2014). Shape, topography, gravity
        anomalies and tidal deformation of Titan. Icarus, 236, 169–177.
        https://doi.org/10.1016/j.icarus.2014.03.018
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10806604/Titan_shape_Mitri2014_unnorm.sh",  # noqa: E501
        known_hash="sha256:19d2dc5be97eff1a60513cfe69dd6ef684ba6e877db60bc0dd5515fa82b13dd5",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )

    return _SHCoeffs.from_file(fname, lmax=lmax, name='Mitri2014_shape',
                               units='m', errors=True, format='dov',
                               encoding='utf-8', normalization='unnorm',
                               csphase=-1)


__all__ = ['Corlies2017_shape', 'Mitri2014_shape']
