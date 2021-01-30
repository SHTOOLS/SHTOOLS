'''
Earth2012: Shape and topography (with respect to mean sea level) of Earth
expanded to degree and order 2160.

shape_air       :  Earth's shape (with water)
shape_bathy     :  Earth's shape (without water)
shape_bathy_bed :  Earth's shape (without water and ice)
shape_ret       :  Earth's rock-equivalent topography as shape model

topo_air        :  Earth's surface (with water)
topo_bathy      :  Earth's surface (without water)
topo_bathy_bed  :  Earth's surface (without water and ice)
ret             :  Earth's rock-equivalent topography

Reference
---------
Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
    isostatic evaluation of new-generation GOCE gravity field models, Journal
    of Geophysical Research: Solid Earth, B05407, doi:10.1029/2011JB008878.
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from ....shclasses import SHCoeffs as _SHCoeffs


def shape_air(lmax=2160):
    '''
    Earth's shape (with water): Harmonic shape model of the interface between
    Earth and its atmosphere, providing radii of the terrain and ice.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.shape_air.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:07b948727c96022f40375d82cb3e505732fb3e1bf72f41b1cc072445dafab059",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='Earth2020.shape_air',
                               units='m', encoding='utf-8')


def shape_bathy(lmax=2160):
    '''
    Earth's shape (without water): Harmonic shape model of the Earth without
    ocean water masses. This model provides radii of the terrain and ice.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.shape_bathy.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:ee03c7addf13f60efd3c928f166f83f5a9f17991c9dd74a88c6c8d9ede8bb15e",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='Earth2012.shape_bathy',
                               units='m', encoding='utf-8')


def shape_bathy_bed(lmax=2160):
    '''
    Earth's shape (without water and ice): Harmonic shape model of the Earth
    without icesheets and without ocean water masses. This model provides radii
    of the terrain over land, of the sea bed over the oceans and inland lakes
    and of bedrock heights over Antarctica and Greenland.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.shape_bathy_bed.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:66ffd58246566b4f62cc1f71145388057a4c2a1a8142f55e089e26a0b2a22d57",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax,
                               name='Earth2012.shape_bathy_bed', units='m',
                               encoding='utf-8')


def shape_ret(lmax=2160):
    '''
    Earth's rock-equivalent topography as shape model: Harmonic shape model of
    Earth's rock-equivalent topography.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.shape_RET2012.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:16c769df9d2c790fb62d1a1e73fd6070c5b5b663c2bd7f68fd21ad4aa8c96a2b",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='Earth2012.shape_ret',
                               units='m', encoding='utf-8')


def topo_air(lmax=2160):
    '''
    Earth's surface (with water): Harmonic model of the interface between
    Earth and its atmosphere, providing heights above mean sea level of the
    terrain and ice over land and zero heights over the oceans.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.topo_air.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:7b68053ba74246f1a755fcce05266a58ab96529b1e48309b98a0e9ba49b4ba3f",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='Earth2012.topo_air',
                               units='m', encoding='utf-8')


def topo_bathy(lmax=2160):
    '''
    Earth's surface (without water): Harmonic model of the Earth's topography
    without ocean water masses. This model provides heights of the terrain and
    of ice over land and bathymetric depths over the oceans, Caspian Sea and
    major inland lakes (Superior, Michigan, Huron, Erie, Ontario and Baikal).

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.topo_bathy.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:997a16f85dafbfd1a0ee2339f503470acea6c8cb8eecdf5f2e057904a97f9718",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='Earth2012.topo_bathy',
                               units='m', encoding='utf-8')


def topo_bathy_bed(lmax=2160):
    '''
    Earth's surface (without water and ice): Harmonic model of the Earth's
    topography without ice sheets and without ocean water masses. This model
    provides heights of the terrain over land, bathymetric depths over the
    oceans, Caspian Sea and major inland lakes, and bedrock heights over
    Antarctica and Greenland.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.topo_bathy_bed.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:f8b0535a76de11de767d0ebb6c032f5983ac48ad29d1750b0d22e15a627f88e1",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax,
                               name='Earth2012.topo_bathy_bed', units='m',
                               encoding='utf-8')


def ret(lmax=2160):
    '''
    Earth's rock-equivalent topography: Harmonic model of Earth's
    rock-equivalent topography.

    Parameters
    ----------
    lmax : int, optional
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Hirt, C., Kuhn, M., Featherstone, W.E., Goettl, F. (2012). Topographic/
        isostatic evaluation of new-generation GOCE gravity field models,
        Journal of Geophysical Research: Solid Earth, B05407,
        doi:10.1029/2011JB008878.
    '''
    fname = _retrieve(
        url="http://ddfe.curtin.edu.au/gravitymodels/Earth2012/topo_shape_to2160/Earth2012.RET2012.SHCto2160.zip",  # noqa: E501
        known_hash="sha256:36b3204d86fa01fa9e8f693e2df8e91905b67619a0192cc0c269be21a7ba5799",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='Earth2012.ret',
                               units='m', encoding='utf-8')


__all__ = ['shape_air', 'shape_bathy', 'shape_bathy_bed', 'shape_ret',
           'topo_air', 'topo_bathy', 'topo_bathy_bed', 'ret']
