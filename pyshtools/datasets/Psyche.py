'''
Datasets related to the asteroid (16) Psyche.

Shape
-----
Shepard2021_shape    :  Wieczorek (2024), Shepard et al. (2021)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHCoeffs as _SHCoeffs


def Shepard2021_shape(lmax=29):
    '''
    Shepard2021 is a spherical harmonic model of the shape of asteroid (16)
    Psyche to degree and order 29. This model was constructed using a least
    squares inversion with the vertice coordinates from the shape model of
    Shepard et al. 2021. The original model makes use of data from Arecibo
    S-band delay-Doppler imaging, Atacama Large Millimeter Array (ALMA)
    imaging, adaptive optics images from Keck and the Very Large Telescope, and
    stellar occultations.

    Parameters
    ----------
    lmax : int, optional, default = 29
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
        (16) Psyche (1.0.0) [Data set]. Zenodo.
        https://doi.org/10.5281/zenodo.10985759
    '''
    if lmax < 0:
        lmax = 29

    fname = _retrieve(
        url="doi:10.5281/zenodo.10985759/Psyche-Shepard2021.sh",
        known_hash="sha256:27107c14e57836d0cde8f11f674f5708fc9aca3449092afd0ca0440a731578fc",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, errors=False,
                               name='Shepard2021_shape (Psyche)',
                               encoding='utf-8'
                               )

__all__ = ['Shepard2021_shape']
