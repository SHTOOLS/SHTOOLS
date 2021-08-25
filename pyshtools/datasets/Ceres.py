'''
Datasets related to the dwarf planet (1) Ceres.

Gravity
-------
CERES18D :  Konopliv et al. (2018)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs


def CERES18D(lmax=18):
    '''
    CERES18D is a JPL 18 degree and order spherical harmonic model of the
    gravitational potential of (1) Ceres.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Park, R.S., Vaughan, A.T., Bills, B.G., Asmar, S.W.,
        Ermakov, A.I., Rambaux, N., Raymond, C.A., Castillo-Rogez, J.C.,
        Russell, C.T., Smith, D.E., Zuber, M.T. (2018). The Ceres gravity
        field, spin pole, rotation period and orbit from the Dawn radiometric
        tracking and optical data, Icarus, 299, 411-429,
        doi:10.1016/j.icarus.2017.08.005.
    '''
    fname = _retrieve(
        url="https://sbnarchive.psi.edu/pds3/dawn/grav/DWNCGRS_2_v3_181005/DATA/SHADR/JGDWN_CER18D_SHA.TAB",  # noqa: E501
        known_hash="sha256:e7ccb1f0c689f77fe5dae4e0bb5d514db1cf5acb5be927bbcaa8576aca153981",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='CERES18D', encoding='utf-8')


__all__ = ['CERES18D']
