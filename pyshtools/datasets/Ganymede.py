'''
Datasets related to the satellite Ganymede.

Gravity
-------
GANYMEDE2022 :  Gomez Casajus et al. (2022)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs


def GANYMEDE2022(lmax=5):
    '''
    GANYMEDE2022 is a JPL spherical harmonic model of the gravitational
    potential of Ganymede to degree and order 5. This model is based on Galileo
    and JUNO extended mission data.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Gomez Casajus, L., Ermakov, A. I., Zannoni, M., Keane, J. T., Stevenson,
        D., Buccino, D. R., Durante, D., Parisi, M., Park, R. S., Tortora,
        P., Bolton, S. J. (2022). Gravity Field of Ganymede After the Juno
        Extended Mission. Geophysical Research Letters, 49(24), e2022GL099475,
        doi:10.1029/2022GL099475.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/7665171/files/Ganymede2022.sh.gz",  # noqa: E501
        known_hash="sha256:593b084cf91673a9093fdfdf657016d52246c6b1a9bf297e42cf543f88fb3b97",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='GANYMEDE2022', encoding='utf-8')


__all__ = ['GANYMEDE2022']
