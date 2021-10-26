'''
Historical datasets related to the magnetic field of Mars.

Magnetic field
-------
FSU50 (Mars Global Surveyor)  :  Cain et al. (2003)
FSU90 (Mars Global Surveyor)  :  Cain et al. (2003)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ....shclasses import SHMagCoeffs as _SHMagCoeffs


def FSU50(lmax=50):
    '''
    FSU50 is a 50th degree and order spherical harmonic model of the
    magnetic potential of Mars. This model makes use of MPO nighttime-only
    data collected from the Mars Global Surveyor mission. The coefficients are
    output in units of nT.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    References
    ----------
    Cain, J. C., Ferguson, B. B., Mozzoni, D. (2003), An n = 90 internal
        potential function of the Martian crustal magnetic field, Journal of
        Geophysical Research: Planets, 108 (E2), doi:10.1029/2000JE001487.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/5503849/files/FSU50.sh.gz?download=1",  # noqa: E501
        known_hash="sha256:ee53cf52869c6144fdeefde96f91c25f3d391768446a3c7c7da1b258eff79efc",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHMagCoeffs.from_file(fname, lmax=lmax, header=True,
                                  file_units='nT', name='FSU50', units='nT',
                                  encoding='utf-8')


def FSU90(lmax=90):
    '''
    FSU90 is a 90th degree and order spherical harmonic model of the
    magnetic potential of Mars. This model makes use of all data collected
    from the Mars Global Surveyor mission. The coefficients are output in units
    of nT.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    References
    ----------
    Cain, J. C., Ferguson, B. B., Mozzoni, D. (2003), An n = 90 internal
        potential function of the Martian crustal magnetic field, Journal of
        Geophysical Research: Planets, 108 (E2), doi:10.1029/2000JE001487.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/5503849/files/FSU90.sh.gz?download=1",  # noqa: E501
        known_hash="sha256:a0c8653c01f06c4af24d011f84ac2f3b119ca19bb79a274f7d5c2536122a4689",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHMagCoeffs.from_file(fname, lmax=lmax, header=True,
                                  file_units='nT', name='FSU90', units='nT',
                                  encoding='utf-8')


__all__ = ['FSU50', 'FSU90']
