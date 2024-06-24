'''
Datasets related to the asteroid (16) Psyche.

Shape
-----
Shepard2017_shape    :  Wieczorek (2024), Shepard et al. (2017)
Shepard2021_shape    :  Wieczorek (2024), Shepard et al. (2021)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHCoeffs as _SHCoeffs


def Shepard2017_shape(lmax=29):
    '''
    Shepard2017 is a spherical harmonic model of the shape of asteroid (16)
    Psyche to degree and order 29. This model was constructed using a least
    squares inversion with the vertice coordinates from the shape model of
    Shepard et al. 2017. The original model makes use of data from Arecibo
    S-band radar data and adaptive-optics images.

    Parameters
    ----------
    lmax : int, optional, default = 29
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Shepard, M. K., Richardson, J., Taylor, P. A., Rodriguez-Ford, L. A.,
        Conrad, A., de Pater, I., Adamkovics, M., de Kleer, K., Males, J. R.,
        Morzinski, K. M., Close, L. M., Kaasalainen, M., Viikinkoski, M.,
        Timerson, B., Reddy, V., Magri, C., Nolan, M. C., Howell, E. S.,
        Benner, L. A. M., Giorgini, J. D., Warner, B. D., Harris, A. W. (2017).
        Radar observations and shape model of asteroid 16 Psyche. Icarus, 281,
        388-403. https://doi.org/10.1016/j.icarus.2016.08.011
    Wieczorek, M. (2024). Spherical harmonic models of the shape of asteroid
        (16) Psyche (1.1) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.12522382
    '''
    if lmax < 0:
        lmax = 29

    fname = _retrieve(
        url="doi:10.5281/zenodo.12522382/Psyche-Shepard2017.sh",
        known_hash="sha256:27107c14e57836d0cde8f11f674f5708fc9aca3449092afd0ca0440a731578fc",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, errors=False,
                               name='Shepard2017_shape (Psyche)',
                               encoding='utf-8'
                               )


def Shepard2021_shape(lmax=10):
    '''
    Shepard2021 is a spherical harmonic model of the shape of asteroid (16)
    Psyche to degree and order 10. This model was constructed using a least
    squares inversion with the vertice coordinates from the shape model of
    Shepard et al. 2021. The original model makes use of data from Arecibo
    S-band delay-Doppler imaging, Atacama Large Millimeter Array (ALMA)
    imaging, adaptive optics images from Keck and the Very Large Telescope, and
    stellar occultations.

    Parameters
    ----------
    lmax : int, optional, default = 10
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Shepard, M. K., de Kleer, K., Cambioni, S., Taylor, P. A., Virkki, A. K.,
        Rívera-Valentin, E. G., Rodriguez Sanchez-Vahamonde, C., Fernanda
        Zambrano-Marin, L., Magri, C., Dunham, D., Moore, J., & Camarca, M.
        (2021). Asteroid 16 Psyche: Shape, Features, and Global Map. The
        Planetary Science Journal, 2(4), 125.
        https://doi.org/10.3847/PSJ/abfdba
    Wieczorek, M. (2024). Spherical harmonic models of the shape of asteroid
        (16) Psyche (1.1) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.12522382
    '''
    if lmax < 0:
        lmax = 10

    fname = _retrieve(
        url="doi:10.5281/zenodo.12522382/Psyche-Shepard2021.sh",
        known_hash="sha256:7d4db7e2f71bb3c7b6d4a881fee532b56bac6204aee8c245884b452d5f27c6d4",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, errors=False,
                               name='Shepard2021_shape (Psyche)',
                               encoding='utf-8'
                               )


__all__ = ['Shepard2017_shape', 'Shepard2021_shape']
