'''
Historical datasets related to the gravity field of Earth's Moon.

Gravity
-------
GLGM2 (Clementine)         :  Lemoine et al. (1997)
GLGM3 (Lunar Prospector)   :  Mazarico et al. (2010)
LP165P (Lunar Prospector)  :  Konopliv et al. (2001)
SGM100i (Kaguya)           :  Goossens et al. (2011)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ....shclasses import SHGravCoeffs as _SHGravCoeffs
from ....constants.Moon import omega as _omega


def GLGM2(lmax=70):
    '''
    GLGM2 is a GSFC 70 degree and order spherical harmonic model of the
    gravitational potential of the Moon derived from Clementine and prior
    mission data. This model applies a Kaula constraint for all degrees.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Lemoine, F., Smith, D., Zuber, M., Neumann, G., Rowlands, D. (1997). A 70th
        degree lunar gravity model (GLGM-2) from Clementine and other tracking
        data. Journal of Geophysical Research 102 (E7), 16,339-16,359,
        doi:10.1029/97JE01418.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/lunar/clem1-gravity-topo-v1/cl_8xxx/gravity/glgm2sh.tab",  # noqa: E501
        known_hash="sha256:94a443391897affa44bdb0c30d5cd5fbf051a6b17df20cff98dc7c8e391734f6",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, header=False,
                                   gm=4902.80295e9, r0=1738.e3,
                                   errors=True, omega=_omega.value,
                                   name='GLGM2', encoding='utf-8')


def GLGM3(lmax=150):
    '''
    GLGM3 is a GSFC 150 degree and order spherical harmonic model of the
    gravitational potential of the Moon derived from Lunar Prospector and prior
    mission data. This model applies a Kaula constraint for all degrees.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Mazarico, E., F. G. Lemoine, S.-C. Han, and D. E. Smith (2010). GLGM-3: A
        degree-150 lunar gravity model from the historical tracking data of
        NASA Moon orbiters, J. Geophys. Res., 115, E05001,
        doi:10.1029/2009JE003472.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/lunar/lp-l-rss-5-glgm3_gravity-v1/lp_1201/data/shadr/gglp_glgm3150_sha.tab",  # noqa: E501
        known_hash="sha256:c8329d8cdfc62e92b9731e0053fe328959bba63444648283bc7df21d42266de4",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, errors=True,
                                   omega=_omega.value, name='GLGM3',
                                   encoding='utf-8')


def LP165P(lmax=165):
    '''
    LP165P is a JPL 165 degree and order spherical harmonic model of the
    gravitational potential of the Moon derived from Lunar Prospector and prior
    mission data. This model applies a Kaula constraint for all degrees.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Konopliv, A. S., S. W. Asmar, E. Carranza, W. L. Sjogren, and D. N. Yuan
        (2001). Recent gravity models as a result of the Lunar Prospector
        mission, Icarus, 150, 1-18, doi:10.1006/icar.2000.6573.
    '''
    fname = _retrieve(
        url="https://pds-geosciences.wustl.edu/lunar/lp-l-rss-5-gravity-v1/lp_1001/sha/jgl165p1.sha",  # noqa: E501
        known_hash="sha256:b75a40be56f9b58bd00536a65c7773db3b6ed17de78fba26b945837fd0d049ee",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, errors=True,
                                   header_units='km', omega=_omega.value,
                                   name='LP165P', encoding='utf-8')


def SGM100I(lmax=100):
    '''
    SGM100I is a JAXA 100 degree and order spherical harmonic model of the
    gravitational potential of the Moon derived from Kaguya and prior
    mission data. This model applies a Kaula constraint for all degrees.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Goossens, S., K. Matsumoto, Q. Liu, F. Kikuchi, K. Sato, H. Hanada,
        Y. Ishihara, H. Noda, N. Kawano, N. Namiki, T. Iwata, F. G. Lemoine,
        D. D. Rowlands, Y. Harada, M. Chen (2011). Lunar gravity field
        determination using SELENE same-beam differential VLBI tracking data.
        J. Geodesy, 85, 205-228, doi:10.1007/s00190-010-0430-2.
    '''
    fname = _retrieve(
        url="https://zenodo.org/record/5233378/files/SGM100i.sh.gz",  # noqa: E501
        known_hash="sha256:0a6e4d8cf3e26f60d61ac194bb63506a4eb91d7dc119a2038cba0acf58849944",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, errors=True,
                                   header_units='m', omega=_omega.value,
                                   name='SGM100I', encoding='utf-8')


__all__ = ['GLGM2', 'GLGM3', 'LP165P', 'SGM100I']
