'''
Datasets related to the satellite Ganymede.

Gravity
-------
Ganymede2022    :  Gomez Casajus et al. (2022)
Anderson1996_1  :  Anderson et al. (1996), encounter 1
Anderson1996_2  :  Anderson et al. (1996), encounter 2

'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..constants.Ganymede import angular_velocity as _omega


def Ganymede2022(lmax=5):
    '''
    Ganymede2022 is a JPL spherical harmonic model of the gravitational
    potential of Ganymede to degree and order 5. This model is based on Galileo
    and JUNO extended mission data.

    Parameters
    ----------
    lmax : int, optional, default = 5
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
        url="doi:10.5281/zenodo.7665171/Ganymede2022.sh.gz",
        known_hash="sha256:593b084cf91673a9093fdfdf657016d52246c6b1a9bf297e42cf543f88fb3b97",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='m',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Ganymede2022', encoding='utf-8',
                                   omega=_omega.value)


def Anderson1996_1(lmax=2):
    '''
    Anderson1996_1 is a JPL spherical harmonic model of the gravitational
    potential of Ganymede to degree and order 2. This model is based on the
    first encounter of Galileo with Ganymede.

    Parameters
    ----------
    lmax : int, optional, default = 2
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Anderson, J. D., Lau, E. L., Sjogren, W. L., Schubert, G., & Moore, W. B.
        (1996). Gravitational constraints on the internal structure of
        Ganymede. Nature, 384(6609), https://doi.org/10.1038/384541a0
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10817282/Anderson1996_Ganymede_1_gravity.sh",
        known_hash="sha256:33e37de6891655c0e57ed7d17fd749525a42cb3c1726078f9b01add2f8c09ceb",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Anderson1996_1 (Ganymede)',
                                   encoding='utf-8', omega=_omega.value,
                                   normalization='unnorm')


def Anderson1996_2(lmax=2):
    '''
    Anderson1996_2 is a JPL spherical harmonic model of the gravitational
    potential of Ganymede to degree and order 2. This model is based on the
    second encounter of Galileo with Ganymede.

    Parameters
    ----------
    lmax : int, optional, default = 2
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Anderson, J. D., Lau, E. L., Sjogren, W. L., Schubert, G., & Moore, W. B.
        (1996). Gravitational constraints on the internal structure of
        Ganymede. Nature, 384(6609), https://doi.org/10.1038/384541a0
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10817282/Anderson1996_Ganymede_2_gravity.sh",
        known_hash="sha256:d2b3f0baa3f84491573dfb520d1b3a2e9e521ead249e54fbf6e280116049495a",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='Anderson1996_2 (Ganymede)',
                                   encoding='utf-8', omega=_omega.value,
                                   normalization='unnorm')


__all__ = ['Anderson1996_1', 'Anderson1996_2', 'Ganymede2022']
