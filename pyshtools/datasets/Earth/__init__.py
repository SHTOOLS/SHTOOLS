'''
Datasets related to the planet Earth.

Topography
----------
Earth2012                  :  Hirt et al. (2012)
Earth2014                  :  Hirt and Rexer (2015)

Gravity
-------
EGM2008                    :  Pavlis et al. (2012); degree 2190
EIGEN_6C4                  :  Förste et al. (2014); degree 2190
GGM05C                     :  Ries et al. (2016); degree 360
GOCO06S                    :  Kvas et al. (2019); degree 300, satellite only
EIGEN_GRGS_RL04_MEAN_FIELD :  Lemoine et al. (2019); degree 300, satellite only
XGM2019E                   :  Zingerle et al. (2019); degree 2190

Magnetic field
--------------
IGRF_13                    :  Thébault et al. (2015)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ...shclasses import SHGravCoeffs as _SHGravCoeffs
from ...shclasses import SHMagCoeffs as _SHMagCoeffs
from ...constant.Earth import omega_egm2008 as _omega

from . import Earth2012
from . import Earth2014


def EGM2008(lmax=2190):
    '''
    EGM2008 is a degree 2190 model of the Earth's gravity field in a
    tide-free system. This model is based on data from altimetry, ground-based
    measurements, and the satellite GRACE. The error coefficients are
    calibrated.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Pavlis, N. K., Holmes, S. A., Kenyon, S. C., Factor, J. K. (2012). The
        development and evaluation of the Earth Gravitational Model 2008
        (EGM2008), Journal of Geophysical Research: Solid Earth, 117, B04406,
        doi:10.1029/2011JB008916.
    '''
    fname = _retrieve(
        url="http://icgem.gfz-potsdam.de/getmodel/gfc/c50128797a9cb62e936337c890e4425f03f0461d7329b09a8cc8561504465340/EGM2008.gfc",  # noqa: E501
        known_hash="sha256:ab5b524da073e63b5bdceb7ca47a0de07a26dd44a1c5798f39fc98dc80af70fd",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, format='icgem',
                                   errors='calibrated', omega=_omega.value)


def EIGEN_6C4(lmax=2190):
    '''
    EIGEN_6C4 is a degree 2190 model of the Earth's gravity field in a
    tide-free system. This model is based on data from altimetry, ground-based
    measurements, and the satellites GOCE, GRACE, and LAGEOS. The error
    coefficients are formal.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Förste, C., Bruinsma, S. L., Abrikosov, O., Lemoine, J.-M., Marty, J. C.,
        Flechtner, F., Balmino, G., Barthelmes, F., Biancale, R. (2014).
        EIGEN-6C4 The latest combined global gravity field model including GOCE
        data up to degree and order 2190 of GFZ Potsdam and GRGS Toulouse, GFZ
        Data Services, http://doi.org/10.5880/icgem.2015.1.
    '''
    fname = _retrieve(
        url="http://icgem.gfz-potsdam.de/getmodel/gfc/7fd8fe44aa1518cd79ca84300aef4b41ddb2364aef9e82b7cdaabdb60a9053f1/EIGEN-6C4.gfc",  # noqa: E501
        known_hash="sha256:2ac274a66a25fb25bdddcec7074669867939345ad003bbc1ece3967450d48dc1",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, format='icgem',
                                   errors='formal', omega=_omega.value,
                                   encoding='iso-8859-1')


def GGM05C(lmax=360):
    '''
    GGM05C is a degree 360 model of the Earth's gravity field in a zero-tide
    system. This model is based on data from altimetry, ground-based
    measurements, and the satellites GOCE and GRACE. The error coefficients are
    calibrated.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Ries, J., Bettadpur, S., Eanes, R., Kang, Z., Ko, U., McCullough, C.,
        Nagel, P., Pie, N., Poole, S., Richter, T., Save, H., Tapley, B.
        (2016). The Combined Gravity Model GGM05C, GFZ Data Services,
        http://doi.org/10.5880/icgem.2016.002.
    '''
    fname = _retrieve(
        url="http://icgem.gfz-potsdam.de/getmodel/gfc/778a683780a5b0ad3163f4772b97b9075a0a13c389d2bd8ea3f891b64cfa383d/GGM05C.gfc",  # noqa: E501
        known_hash="sha256:efb1781e748969dab2721f45365cd060b5963ffe4fd497006cc786f9e5241bf2",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, format='icgem',
                                   errors='calibrated', omega=_omega.value)


def GOCO06S(lmax=300):
    '''
    GOCO06S is a degree 300 model of the Earth's gravity field in a zero-tide
    system. This model is based solely on satellite data. The error
    coefficients are formal.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Kvas, A., Mayer-Gürr, T., Krauss, S., Brockmann, J. M., Schubert, T.,
        Schuh, W.-D., Pail, R., Gruber, T., Jäggi, A., Meyer, U. (2019).
        The satellite-only gravity field model GOCO06s, GFZ Data Services,
        http://doi.org/10.5880/ICGEM.2019.002.
    '''
    fname = _retrieve(
        url="http://icgem.gfz-potsdam.de/getmodel/gfc/32ec2884630a02670476f752d2a2bf1c395d8c8d6d768090ed95b4f04b0d5863/GOCO06s.gfc",  # noqa: E501
        known_hash="sha256:351d9d20b84cd2c0f52ce77146b1e3b774f408200b579ffaf98593cf3d271819",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, format='icgem',
                                   errors='formal', omega=_omega.value)


def EIGEN_GRGS_RL04_MEAN_FIELD(epoch=None, lmax=300):
    '''
    EIGEN_GRGS_RL04_MEAN_FIELD is a degree 300 model of the Earth's gravity
    field in a tide-free system. This model is based solely on satellite
    data. The error coefficients are formal.

    Parameters
    ----------
    epoch : str or float
        The epoch time to calculate time-variable coefficients in YYYYMMDD.DD
        format.
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Lemoine, J.-M., Bourgogne, S., Biancale, R., Reinquin, F. (2019). The new
        time-variable gravity field model for POD of altimetric satellites
        based on GRACE+SLR RL04 from CNES/GRGS, IDS Workshop, 24-29 September,
        Ponta Delgada, Portugal.
    '''
    fname = _retrieve(
        url="http://icgem.gfz-potsdam.de/getmodel/gfc/96e984f23a31505411d943341767e9ff082ad2c77fbbdcba1be2e393cf35110a/EIGEN-GRGS.RL04.MEAN-FIELD.gfc",  # noqa: E501
        known_hash="sha256:ea5d99d84da42d0eedccb381c282604a1a42fd64d1f9e25bba796b7e8838732f",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, format='icgem',
                                   errors='formal', omega=_omega.value,
                                   epoch=epoch)


def XGM2019E(lmax=2190):
    '''
    XGM2019E is a degree 2190 experimental model of the Earth's gravity field
    in a zero-tide system. This model is based on data from altimetry,
    ground-based measurements, topography, and satellites. The error
    coefficients are formal.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Zingerle, P., Pail, R., Gruber, T., Oikonomidou, X. (2019). The
        experimental gravity field model XGM2019e, GFZ Data Services,
        http://doi.org/10.5880/ICGEM.2019.007.
    '''
    fname = _retrieve(
        url="http://icgem.gfz-potsdam.de/getmodel/gfc/eeb03971cf6e533e6eeb6b010336463286dcda0846684248d5530acf8e800055/XGM2019e_2159.gfc",  # noqa: E501
        known_hash="sha256:ee24f7cf68f45b44cdbf3ccbf57e7382c7667869356061419f00f2ca951861e4",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, format='icgem',
                                   errors='formal', omega=_omega.value)


def IGRF_13(lmax=13, year=2020.):
    '''
    IGRF-13 is a degree 13 time variable model of the Earth's main magnetic
    field that is valid between 1900 and 2020. Coefficients are provided in 5
    year intervals, and for a given year, the values of the coefficients are
    interpolated linearly between adjacent entries. For years between 2020 and
    2025, the coefficients are extrapolated using the provided secular
    variation.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.
    year : float, optional, default = 2020.
        The year to compute the coefficients.

    Reference
    ---------
    Thébault, E., and 47 coauthors (2015). International Geomagnetic Reference
        Field: the 12th generation, Earth, Planets and Space, 67, 79,
        doi:10.1186/s40623-015-0228-9.
    '''
    fname = _retrieve(
        url="https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt",
        known_hash="sha256:460b8d8beb9b4df84febe4f0b639f0dd54dccfe8ff0970616287b015fa721425",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHMagCoeffs.from_file(fname, format='igrf', r0=6371.2e3,
                                  lmax=lmax, year=year)


__all__ = ['Earth2012', 'Earth2014', 'EGM2008', 'EIGEN_6C4',
           'GGM05C', 'GOCO06S', 'EIGEN_GRGS_RL04_MEAN_FIELD', 'XGM2019E',
           'IGRF_13']
