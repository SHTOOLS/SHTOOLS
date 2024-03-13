'''
Datasets related to Saturn's moon Enceladus.

Shape
-----
JPL_SPC_shape     :  Wieczorek (2024)

Gravity
-------
Iess2014_gravity  :  Iess et al. (2014)
Park2020_gravity  :  Park et al. (2024)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import create as _create
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHCoeffs as _SHCoeffs
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Enceladus import omega as _omega


def JPL_SPC_shape(lmax=719):
    '''
    JPL_SPC_shape is a spherical harmonic model of the shape of Saturn's moon
    Enceladus based on stereo photoclinometric data obtained by the Cassini
    mission. The maximum spherical harmonic degree of the model is 1023, which
    has an effective spatial resolution of 11.3 pixels per degree. One lower
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
        shape of Enceladus [JPL SPC] (1.0.0) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.10813481
    Park, R. S., Mastrodemos, N., Jacobson, R. A., Berne, A., Vaughan, A. T.,
        Hemingway, D. J., Leonard, E. J., Castillo-Rogez, J. C., Cockell,
        C. S., Keane, J. T., Konopliv, A. S., Nimmo, F., Riedel, J. E.,
        Simons, M., & Vance, S. (2024). The Global Shape, Gravity Field, and
        Libration of Enceladus. Journal of Geophysical Research: Planets,
        129(1), e2023JE008054. https://doi.org/10.1029/2023JE008054
    '''
    archive = _create(
        path=_os_cache('pyshtools'),
        base_url="doi:10.5281/zenodo.10813481",
        registry={
            "Enceladus_JPL_SPC_shape_1023.sh.gz": "sha256:8094ca35956f318c644497c90afc2e31135596128f3dd8eeefdfdc085367b7c6",  # noqa: E501
            "Enceladus_JPL_SPC_shape_719.sh.gz": "sha256:dd126110b49c28c4a5df466bcc50257ea2e467973320e7ccbffb3900f5a45e41",  # noqa: E501
            },
        )

    if lmax < 0:
        lmax = 1023

    if lmax >= 0 and lmax <= 719:
        fname = archive.fetch("Enceladus_JPL_SPC_shape_719.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
    else:
        fname = archive.fetch("Enceladus_JPL_SPC_shape_1023.sh.gz",
                              downloader=_DOIDownloader(progressbar=True))
        lmax = min(lmax, 1023)

    return _SHCoeffs.from_file(fname, lmax=lmax,
                               name='JPL_SPC_shape (Enceladus)',
                               units='m', format='bshc')


def Iess2014_gravity(lmax=3):
    '''
    Iess2014_gravity is a degree and order 3 spherical harmonic model of the
    gravitational potential of Saturn's moon Enceladus. This model corresponds
    to SOL1 in Iess et al. (2014).

    Parameters
    ----------
    lmax : int, optional, default = 3
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Iess, L., Stevenson, D. J., Parisi, M., Hemingway, D., Jacobson, R. A.,
        Lunine, J. I., Nimmo, F., Armstrong, J. W., Asmar, S. W., Ducci, M.,
        & Tortora, P. (2014). The Gravity Field and Interior Structure of
        Enceladus. Science, 344(6179), 78–80.
        https://doi.org/10.1126/science.1250551
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10809142/Iess2014.sh",  # noqa: E501
        known_hash="sha256:f393e0ab8a7c0420907ed1f9e33b51c8bf38589e8e67ebdd3a8ab35ba1ad6988",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Iess2014_gravity (Enceladus)',
                                   encoding='utf-8',
                                   omega=_omega.value,
                                   normalization='unnorm')


def Park2024_gravity(lmax=3):
    '''
    Park2014_gravity is a degree and order 3 spherical harmonic model of the
    gravitational potential of Saturn's moon Enceladus. This model corresponds
    to Case 2 in Park et al. (2024).

    Parameters
    ----------
    lmax : int, optional, default = 3
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Park, R. S., Mastrodemos, N., Jacobson, R. A., Berne, A., Vaughan, A. T.,
        Hemingway, D. J., Leonard, E. J., Castillo-Rogez, J. C., Cockell,
        C. S., Keane, J. T., Konopliv, A. S., Nimmo, F., Riedel, J. E.,
        Simons, M., & Vance, S. (2024). The Global Shape, Gravity Field, and
        Libration of Enceladus. Journal of Geophysical Research: Planets,
        129(1), e2023JE008054. https://doi.org/10.1029/2023JE008054
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10809142/Park2024.sh",  # noqa: E501
        known_hash="sha256:f1b29bec0132f008032a8dbf07366a64cd33c404787ce0839ed1799d4863318f",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Park2024_gravity (Enceladus)',
                                   encoding='utf-8',
                                   omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['JPL_SPC_shape', 'Iess2014_gravity', 'Park2024_gravity']
