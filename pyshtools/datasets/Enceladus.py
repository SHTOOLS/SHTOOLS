'''
Datasets related to Saturn's moon Enceladus.

Gravity
-------
Iess2014_gravity  :  Iess et al. (2014)
Park2020_gravity  :  Park et al. (2024)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
# from ..shclasses import SHCoeffs as _SHCoeffs
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Enceladus import omega as _omega


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
                                   name='Iess2014_gravity', encoding='utf-8',
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
                                   name='Park2024_gravity', encoding='utf-8',
                                   omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Iess2014_gravity', 'Park2024_gravity']
