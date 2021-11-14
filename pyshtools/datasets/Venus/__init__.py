'''
Datasets related to the planet Venus.

Topography
----------
VenusTopo719 :  Wieczorek (2015)

Gravity
-------
MGNP180U     :  Konopliv et al. (1999)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ...shclasses import SHCoeffs as _SHCoeffs
from ...shclasses import SHGravCoeffs as _SHGravCoeffs
from ...constants.Venus import omega as _omega
from . import historical  # noqa: F401


def VenusTopo719(lmax=719):
    '''
    VenusTopo719 is a 719 degree and order spherical harmonic model of the
    shape of the planet Venus. The coefficients are in units of meters.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Wieczorek, M.A. (2015). Gravity and Topography of the Terrestrial Planets,
        Treatise on Geophysics, 2nd edition, Oxford, 153-193,
        doi:10.1016/B978-0-444-53802-4.00169-X.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/3870926/files/VenusTopo719.shape.gz",
        known_hash="sha256:9fcb04fb21eb7090df68e42458f6d7495a27ff62687b46534057ed608786cf3b",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='VenusTopo719',
                               units='m', encoding='utf-8')


def MGNP180U(lmax=180):
    '''
    MGNP180U is a JPL 180 degree and order spherical harmonic model of the
    gravitational potential of Venus.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Banerdt, W.B., Sjogren, W.L. (1999). Venus gravity: 180th
        degree and order model, Icarus, 139, 3-18, doi:10.1006/icar.1999.6086.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/shgj180u.a01",  # noqa: E501
        known_hash="sha256:d59e2bb90c104ca1157681454a4c016ed4d1cb0e496861b1d3b35829403cf53b",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, skip=236, header=True,
                                   r0_index=1, gm_index=2, header_units='km',
                                   errors=True, omega=_omega.value,
                                   name='MGNP180U', encoding='utf-8')


__all__ = ['VenusTopo719', 'MGNP180U', 'historical']
