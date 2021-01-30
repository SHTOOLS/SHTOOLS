'''
Datasets related to the asteroid (4) Vesta.

Gravity
-------
VESTA20H :  Konopliv et al. (2014)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs


def VESTA20H(lmax=20):
    '''
    VESTA20H is a JPL 20 degree and order spherical harmonic model of the
    gravitational potential of (4) Vesta.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A.S., Asmar, S.W., Park, R.S., Bills, B.G., Centinello, F.,
        Chamberlin, A.B., Ermakov, A., Gaskell, R.W., Rambaux, N.,
        Raymond, C.A., Russell, C.T., Smith, D.E., Tricarico, P., Zuber, M.T.
        (2014). The Vesta gravity field, spin pole and rotation period,
        landmark positions, and ephemeris from the Dawn tracking and optical
        data. Icarus, 240, 103-117, doi:10.1016/j.icarus.2013.09.005.
    '''
    fname = _retrieve(
        url="https://sbnarchive.psi.edu/pds3/dawn/grav/DWNVGRS_2/DATA/SHADR/JGDWN_VES20H_SHA.TAB",  # noqa: E501
        known_hash="sha256:389bb30882b43b36877e9ff8bf54c55a210d7ba98cab88e9c4a1fa83ad8484b7",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header_units='km',
                                   r0_index=0, gm_index=1, errors=True,
                                   name='VESTA20H', encoding='utf-8')


__all__ = ['VESTA20H']
