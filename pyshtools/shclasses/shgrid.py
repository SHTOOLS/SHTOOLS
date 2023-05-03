"""
    Spherical Harmonic Grid classes
"""
import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable
import copy as _copy
import xarray as _xr
import tempfile as _tempfile
from ..backends import backend_module
from ..backends import preferred_backend
from ..backends import shtools as _shtools

try:
    import cartopy.crs as _ccrs
    from cartopy.mpl.ticker import LongitudeFormatter as _LongitudeFormatter
    from cartopy.mpl.ticker import LatitudeFormatter as _LatitudeFormatter
    _cartopy_module = True
except ModuleNotFoundError:
    _cartopy_module = False

try:
    import pygmt as _pygmt
    _pygmt_module = True
except ModuleNotFoundError:
    _pygmt_module = False


class SHGrid(object):
    """
    Class for spatial gridded data on the sphere.

    Grids can be initialized from:

        x = SHGrid.from_array(array)
        x = SHGrid.from_xarray(data_array)
        x = SHGrid.from_netcdf(netcdf)
        x = SHGrid.from_file('fname.dat')
        x = SHGrid.from_zeros(lmax)
        x = SHGrid.from_cap(theta, clat, clon, lmax)
        x = SHGrid.from_ellipsoid(lmax, a, b, c)

    The class instance defines the following class attributes:

    data       : Gridded array of the data.
    nlat, nlon : The number of latitude and longitude bands in the grid.
    n          : The number of samples in latitude for 'DH' grids.
    lmax       : The maximum spherical harmonic degree that can be resolved
                 by the grid sampling.
    sampling   : The longitudinal sampling for Driscoll and Healy grids. Either
                 1 for equally sampled grids (nlat=nlon) or 2 for equally
                 spaced grids in degrees.
    kind       : Either 'real' or 'complex' for the data type.
    grid       : Either 'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-
                 Legendre Quadrature grids.
    units      : The units of the gridded data.
    zeros      : The cos(colatitude) nodes used with Gauss-Legendre
                 Quadrature grids. Default is None.
    weights    : The latitudinal weights used with Gauss-Legendre
                 Quadrature grids. Default is None.
    extend     : True if the grid contains the redundant column for 360 E and
                 (for 'DH' grids) the unnecessary row for 90 S.

    Each class instance provides the following methods:

    to_array()       : Return the raw gridded data as a numpy array.
    to_xarray()      : Return the gridded data as an xarray DataArray.
    to_file()        : Save gridded data to a text or binary file.
    to_netcdf()      : Return the gridded data as a netcdf formatted file or
                       object.
    to_real()        : Return a new SHGrid class instance of the real component
                       of the data.
    to_imag()        : Return a new SHGrid class instance of the imaginary
                       component of the data.
    lats()           : Return a vector containing the latitudes of each row
                       of the gridded data.
    lons()           : Return a vector containing the longitudes of each column
                       of the gridded data.
    histogram()      : Return an area-weighted histogram of the gridded data.
    expand()         : Expand the grid into spherical harmonics.
    max()            : Return the maximum value of data using numpy.max().
    min()            : Return the minimum value of data using numpy.min().
    copy()           : Return a copy of the class instance.
    plot()           : Plot the data.
    plotgmt()        : Plot projected data using the generic mapping tools
                       (GMT).
    plot3d()         : Plot a 3-dimensional representation of the data.
    plot_histogram() : Plot a histogram of the area-weighted gridded data.
    info()           : Print a summary of the data stored in the SHGrid
                       instance.
    """

    def __init__():
        """Unused constructor of the super class."""
        print('Initialize the class using one of the class methods:\n'
              '>>> pyshtools.SHGrid.from_array\n'
              '>>> pyshtools.SHGrid.from_xarray\n'
              '>>> pyshtools.SHGrid.from_netcdf\n'
              '>>> pyshtools.SHGrid.from_file\n'
              '>>> pyshtools.SHGrid.from_zeros\n'
              '>>> pyshtools.SHGrid.from_cap\n'
              '>>> pyshtools.SHGrid.from_ellipsoid\n')

    # ---- Factory methods ----
    @classmethod
    def from_array(self, array, grid='DH', units=None, copy=True):
        """
        Initialize the class instance from an input array.

        Usage
        -----
        x = SHGrid.from_array(array, [grid, units, copy])

        Returns
        -------
        x : SHGrid class instance

        Parameters
        ----------
        array : ndarray, shape (nlat, nlon)
            2-D numpy array of the gridded data, where nlat and nlon are the
            number of latitudinal and longitudinal bands, respectively.
        grid : str, optional, default = 'DH'
            'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-Legendre
            Quadrature grids, respectively.
        units : str, optional, default = None
            The units of the gridded data.
        copy : bool, optional, default = True
            If True (default), make a copy of array when initializing the class
            instance. If False, initialize the class instance with a reference
            to array.
        """
        if _np.iscomplexobj(array):
            kind = 'complex'
        else:
            kind = 'real'

        if type(grid) != str:
            raise ValueError('grid must be a string. Input type is {:s}.'
                             .format(str(type(grid))))

        if grid.upper() not in set(['DH', 'GLQ']):
            raise ValueError(
                "grid must be 'DH' or 'GLQ'. Input value is {:s}."
                .format(repr(grid))
                )

        for cls in self.__subclasses__():
            if cls.istype(kind) and cls.isgrid(grid):
                return cls(array, units=units, copy=copy)

    @classmethod
    def from_zeros(self, lmax, grid='DH', kind='real', sampling=2,
                   units=None, extend=True, empty=False):
        """
        Initialize the class instance using an array of zeros.

        Usage
        -----
        x = SHGrid.from_zeros(lmax, [grid, kind, sampling, units, extend,
                                     empty])

        Returns
        -------
        x : SHGrid class instance

        Parameters
        ----------
        lmax : int
            The maximum spherical harmonic degree resolvable by the grid.
        grid : str, optional, default = 'DH'
            'DH' or 'GLQ' for Driscoll and Healy grids or Gauss Legendre
            Quadrature grids, respectively.
        kind : str, optional, default = 'real'
            Either 'real' or 'complex' for the data type.
        sampling : int, optional, default = 2
            The longitudinal sampling for Driscoll and Healy grids. Either 1
            for equally sampled grids (nlong=nlat) or 2 for equally spaced
            grids in degrees (nlong=2*nlat with extend=False or nlong=2*nlat-1
            with extend=True).
        units : str, optional, default = None
            The units of the gridded data.
        extend : bool, optional, default = True
            If True, include the longitudinal band for 360 E (DH and GLQ grids)
            and the latitudinal band for 90 S (DH grids only).
        empty : bool, optional, default = False
            If True, create the data array using numpy.empty() and do not
            initialize with zeros.
        """
        if type(grid) != str:
            raise ValueError('grid must be a string. Input type is {:s}.'
                             .format(str(type(grid))))

        if grid.upper() not in set(['DH', 'GLQ']):
            raise ValueError("grid must be 'DH' or 'GLQ'. " +
                             "Input value is {:s}.".format(repr(grid)))

        if grid.upper() == 'DH':
            nlat = 2 * lmax + 2
            if sampling == 1:
                nlon = nlat
            else:
                nlon = nlat * 2
            if extend:
                nlat += 1
                nlon += 1
        elif grid.upper() == 'GLQ':
            nlat = lmax + 1
            nlon = 2 * nlat - 1
            if extend:
                nlon += 1

        if kind == 'real':
            if empty:
                array = _np.empty((nlat, nlon), dtype=_np.float64)
            else:
                array = _np.zeros((nlat, nlon), dtype=_np.float64)
        else:
            if empty:
                array = _np.empty((nlat, nlon), dtype=_np.complex128)
            else:
                array = _np.zeros((nlat, nlon), dtype=_np.complex128)

        for cls in self.__subclasses__():
            if cls.istype(kind) and cls.isgrid(grid):
                return cls(array, units=units, copy=False)

    @classmethod
    def from_ellipsoid(self, lmax, a, b=None, c=None, grid='DH', kind='real',
                       sampling=2, units=None, extend=True):
        """
        Initialize the class instance with a triaxial ellipsoid whose principal
        axes are aligned with the x, y, and z axes.

        Usage
        -----
        x = SHGrid.from_ellipsoid(lmax, a, [b, c, grid, kind, sampling,
                                            units, extend])

        Returns
        -------
        x : SHGrid class instance

        Parameters
        ----------
        a : float
            Length of the principal axis aligned with the x axis.
        b : float, optional, default = a
            Length of the principal axis aligned with the y axis.
        c : float, optional, default = b
            Length of the principal axis aligned with the z axis.
        lmax : int
            The maximum spherical harmonic degree resolvable by the grid.
        grid : str, optional, default = 'DH'
            'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-Legendre
            Quadrature grids, respectively.
        kind : str, optional, default = 'real'
            Either 'real' or 'complex' for the data type.
        sampling : int, optional, default = 2
            The longitudinal sampling for Driscoll and Healy grids. Either 1
            for equally sampled grids (nlong=nlat) or 2 for equally spaced
            grids in degrees (nlong=2*nlat with extend=False or nlong=2*nlat-1
            with extend=True).
        units : str, optional, default = None
            The units of the gridded data.
        extend : bool, optional, default = True
            If True, include the longitudinal band for 360 E (DH and GLQ grids)
            and the latitudinal band for 90 S (DH grids only).
        """
        temp = self.from_zeros(lmax, grid=grid, kind=kind, sampling=sampling,
                               units=units, extend=extend, empty=True)
        if c is None and b is None:
            temp.data[:, :] = a
        elif c is not None and b is None:
            for ilat, lat in enumerate(temp.lats()):
                temp.data[ilat, :] = 1. / _np.sqrt(
                    _np.cos(_np.deg2rad(lat))**2 / a**2 +
                    _np.sin(_np.deg2rad(lat))**2 / c**2
                    )
        else:
            if c is None:
                c = b
            cos2 = _np.cos(_np.deg2rad(temp.lons()))**2
            sin2 = _np.sin(_np.deg2rad(temp.lons()))**2
            for ilat, lat in enumerate(temp.lats()):
                temp.data[ilat, :] = 1. / _np.sqrt(
                    _np.cos(_np.deg2rad(lat))**2 * cos2 / a**2 +
                    _np.cos(_np.deg2rad(lat))**2 * sin2 / b**2 +
                    _np.sin(_np.deg2rad(lat))**2 / c**2
                    )

        return temp

    @classmethod
    def from_cap(self, theta, clat, clon, lmax, grid='DH', kind='real',
                 sampling=2, degrees=True, units=None, extend=True):
        """
        Initialize the class instance with an array equal to unity within
        a spherical cap and zero elsewhere.

        Usage
        -----
        x = SHGrid.from_cap(theta, clat, clon, lmax, [grid, kind, sampling,
                            degrees, units, extend])

        Returns
        -------
        x : SHGrid class instance

        Parameters
        ----------
        theta : float
            The angular radius of the spherical cap, default in degrees.
        clat, clon : float
            Latitude and longitude of the center of the rotated spherical cap
            (default in degrees).
        lmax : int
            The maximum spherical harmonic degree resolvable by the grid.
        grid : str, optional, default = 'DH'
            'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-Legendre
            Quadrature grids, respectively.
        kind : str, optional, default = 'real'
            Either 'real' or 'complex' for the data type.
        sampling : int, optional, default = 2
            The longitudinal sampling for Driscoll and Healy grids. Either 1
            for equally sampled grids (nlong=nlat) or 2 for equally spaced
            grids in degrees (nlong=2*nlat with extend=False or nlong=2*nlat-1
            with extend=True).
        degrees : bool, optional = True
            If True, theta, clat, and clon are in degrees.
        units : str, optional, default = None
            The units of the gridded data.
        extend : bool, optional, default = True
            If True, include the longitudinal band for 360 E (DH and GLQ grids)
            and the latitudinal band for 90 S (DH grids only).
        """
        temp = self.from_zeros(lmax, grid=grid, kind=kind, sampling=sampling,
                               units=units, extend=extend)

        if degrees is True:
            theta = _np.deg2rad(theta)
            clat = _np.deg2rad(clat)
            clon = _np.deg2rad(clon)

        # Set array equal to 1 within the cap
        lats = temp.lats(degrees=False)
        lons = temp.lons(degrees=False)
        imin = _np.inf
        imax = 0
        for i, lat in enumerate(lats):
            if lat <= clat + theta:
                if i <= imin:
                    imin = i
            if lat >= clat - theta:
                if i >= imax:
                    imax = i

        x = _np.cos(clat) * _np.cos(clon)
        y = _np.cos(clat) * _np.sin(clon)
        z = _np.sin(clat)

        coslon = _np.cos(lons)
        sinlon = _np.sin(lons)
        costheta = _np.cos(theta)
        for i in range(imin, imax+1):
            coslat = _np.cos(lats[i])
            sinlat = _np.sin(lats[i])
            for j in range(0, temp.nlon):
                dist = coslat * (x * coslon[j] + y * sinlon[j]) + z * sinlat
                if dist >= costheta:
                    # ie. _np.arccos(dist) <= theta
                    # since 0 <= theta <= pi/2 and 0 <= dist <= 1
                    # cos is decreasing
                    temp.data[i, j] = 1.

        return temp

    @classmethod
    def from_file(self, fname, binary=False, grid='DH', units=None, **kwargs):
        """
        Initialize the class instance from gridded data in a file.

        Usage
        -----
        x = SHGrid.from_file(fname, [binary, grid, units, **kwargs])

        Returns
        -------
        x : SHGrid class instance

        Parameters
        ----------
        fname : str
            The filename containing the gridded data. For text files (default)
            the file is read using the numpy routine loadtxt(), whereas for
            binary files, the file is read using numpy.load(). For Driscoll and
            Healy grids, the dimensions of the array must be nlon=nlat,
            nlon=2*nlat or nlon=2*nlat-1. For Gauss-Legendre Quadrature grids,
            the dimensions of the array must be nlon=2*nlat-1 or nlon=2*nlat.
            For text files, if the filename ends in '.gz', the file will be
            decompressed using gzip.
        binary : bool, optional, default = False
            If False, read a text file. If True, read a binary 'npy' file.
        grid : str, optional, default = 'DH'
            'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-Legendre
            Quadrature grids, respectively.
        units : str, optional, default = None
            The units of the gridded data.
        **kwargs : keyword arguments, optional
            Keyword arguments of numpy.loadtxt() or numpy.load().
        """
        if binary is False:
            data = _np.loadtxt(fname, **kwargs)
        elif binary is True:
            data = _np.load(fname, **kwargs)
        else:
            raise ValueError('binary must be True or False. '
                             'Input value is {:s}.'.format(binary))

        return self.from_array(data, grid=grid, units=units, copy=False)

    @classmethod
    def from_xarray(self, data_array, grid='DH', units=None):
        """
        Initialize the class instance from an xarray DataArray object.

        Usage
        -----
        x = SHGrid.from_xarray(data_array, [grid])

        Returns
        -------
        x : SHGrid class instance

        Parameters
        ----------
        xarray : xarray DataArray
            The xarray DataArray containing the gridded data. For Driscoll and
            Healy grids, the dimensions of the array must be nlon=nlat,
            nlon=2*nlat or nlon=2*nlat-1. For Gauss-Legendre Quadrature grids,
            the dimensions of the array must be nlon=2*nlat-1 or nlon=2*nlat.
        grid : str, optional, default = 'DH'
            'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-Legendre
            Quadrature grids, respectively.
        units : str, optional, default = None
            The units of the gridded data.
        """
        try:
            units = data_array.units
        except:
            pass

        return self.from_array(data_array.values, grid=grid, units=units)

    @classmethod
    def from_netcdf(self, netcdf, grid='DH', units=None):
        """
        Initialize the class instance from a netcdf formatted file or object.

        Usage
        -----
        x = SHGrid.from_netcdf(netcdf, [grid])

        Returns
        -------
        x : SHGrid class instance

        Parameters
        ----------
        netcdf : str or netcdf object
            The name of a netcdf file or object.
        grid : str, optional, default = 'DH'
            'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-Legendre
            Quadrature grids, respectively.
        units : str, optional, default = None
            The units of the gridded data.
        """
        data_array = _xr.open_dataarray(netcdf)

        try:
            units = data_array.units
        except:
            pass

        return self.from_array(data_array.values, grid=grid, units=units)

    # ---- I/O methods ----
    def copy(self):
        """
        Return a deep copy of the class instance.

        Usage
        -----
        copy = x.copy()
        """
        return _copy.deepcopy(self)

    def to_file(self, filename, binary=False, **kwargs):
        """
        Save gridded data to a file.

        Usage
        -----
        x.to_file(filename, [binary, **kwargs])

        Parameters
        ----------
        filename : str
            Name of output file. For text files (default), the file will be
            compressed using gzip if filename ends in '.gz.'
        binary : bool, optional, default = False
            If False, save as text using numpy.savetxt(). If True, save as a
            'npy' binary file using numpy.save().
        **kwargs : keyword arguments, optional
            Keyword arguments of numpy.savetxt() and numpy.save().
        """
        if binary is False:
            _np.savetxt(filename, self.data, **kwargs)
        elif binary is True:
            _np.save(filename, self.data, **kwargs)
        else:
            raise ValueError('binary must be True or False. '
                             'Input value is {:s}.'.format(binary))

    def to_xarray(self, title=None, comment='pyshtools grid',
                  long_name=None, units=None):
        """
        Return the gridded data as an xarray DataArray.

        Usage
        -----
        x.to_xarray([title, comment, long_name, units])

        Parameters
        ----------
        title : str, optional, default = None
            Title of the dataset.
        comment : str, optional, default = 'pyshtools grid'
            Additional information about how the data were generated.
        long_name : str, optional, default = None
            A long descriptive name of the gridded data, used to label a
            colorbar.
        units : str, optional, default = None
            Units of the gridded data, used to label a colorbar.
        """
        attrs = {'actual_range': [self.min(), self.max()],
                 'comment': comment,
                 'nlat': self.nlat,
                 'nlon': self.nlon,
                 'lmax': self.lmax,
                 'kind': self.kind,
                 'grid': self.grid,
                 'extend': repr(self.extend)
                 }
        if self.grid == 'GLQ':
            attrs['zeros'] = self.zeros
            attrs['weights'] = self.weights
        else:
            attrs['sampling'] = self.sampling
        if title is not None:
            attrs['title'] = title
        if long_name is not None:
            attrs['long_name'] = long_name
        if units is None:
            units = self.units
        if units is not None:
            attrs['units'] = units

        da = _xr.DataArray(self.to_array(),
                           coords=[('lat', self.lats(),
                                    {'long_name': 'latitude',
                                     'units': 'degrees_north',
                                     'actual_range': [self.lats()[0],
                                                      self.lats()[-1]]}),
                                   ('lon', self.lons(),
                                    {'long_name': 'longitude',
                                     'units': 'degrees_east',
                                     'actual_range': [self.lons()[0],
                                                      self.lons()[-1]]})],
                           attrs=attrs)
        if _pygmt_module:
            da.gmt.registration = 0
            da.gmt.gtype = 1
        return da

    def to_netcdf(self, filename=None, title=None, description=None,
                  comment='pyshtools grid', name='data',
                  long_name=None, units=None, dtype='d'):
        """
        Return the gridded data as a netcdf formatted file or object.

        Usage
        -----
        x.to_netcdf([filename, title, description, comment, name, long_name,
                     units, dtype])

        Parameters
        ----------
        filename : str, optional, default = None
            Name of output file.
        title : str, optional, default = None
            Title of the dataset.
        description : str, optional, default = None
            Description of the dataset ('Remark' in gmt grd files).
        comment : str, optional, default = 'pyshtools grid'
            Additional information about how the data were generated.
        name : str, optional, default = 'data'
            Name of the data array.
        long_name : str, optional, default = None
            A long descriptive name of the gridded data.
        units : str, optional, default = None
            Units of the gridded data.
        dtype : str, optional, default = 'd'
            Data type of the output array. Either 'f' or 'd' for single or
            double precision floating point, respectively.
        """
        if self.kind == 'complex':
            raise RuntimeError('netcdf files do not support complex data '
                               'formats.')

        _data = self.to_xarray(title=title, comment=comment,
                               long_name=long_name, units=units)

        if dtype == 'f':
            _data.values = _data.values.astype(_np.float32)
        elif dtype != 'd':
            raise ValueError("dtype must be either 'f' or 'd' for single or "
                             "double precision floating point.")

        attrs = {}
        if title is not None:
            attrs['title'] = title
        if description is not None:
            attrs['description'] = description
        if comment is not None:
            attrs['comment'] = comment
        if units is None:
            units = self.units
        if units is not None:
            attrs['units'] = units

        _dataset = _xr.Dataset({name: _data}, attrs=attrs)

        if filename is None:
            return _dataset.to_netcdf()
        else:
            _dataset.to_netcdf(filename)

    def to_array(self):
        """
        Return the raw gridded data as a numpy array.

        Usage
        -----
        grid = x.to_array()

        Returns
        -------
        grid : ndarray, shape (nlat, nlon)
            2-D numpy array of the gridded data.
        """
        return _np.copy(self.data)

    def to_real(self):
        """
        Return a new SHGrid class instance of the real component of the data.

        Usage
        -----
        grid = x.to_real()

        Returns
        -------
        grid : SHGrid class instance
        """
        return SHGrid.from_array(self.to_array().real, grid=self.grid,
                                 units=self.units, copy=False)

    def to_imag(self):
        """
        Return a new SHGrid class instance of the imaginary component of the
        data.

        Usage
        -----
        grid = x.to_imag()

        Returns
        -------
        grid : SHGrid class instance
        """
        return SHGrid.from_array(self.to_array().imag, grid=self.grid,
                                 units=self.units, copy=False)

    def info(self):
        """
        Print a summary of the data stored in the SHGrid instance.

        Usage
        -----
        x.info()
        """
        print(repr(self))

    # -------------------------------------------------------------------------
    #    Mathematical operators
    #
    #    All operations ignore the units of the coefficients, with the
    #    exception of multiplying and dividing by a scalar.
    # -------------------------------------------------------------------------
    def min(self):
        """
        Return the minimum value of self.data using numpy.min().

        Usage
        -----
        x.min()
        """
        return _np.min(self.data)

    def max(self):
        """
        Return the maximum value of self.data using numpy.max().

        Usage
        -----
        x.max()
        """
        return _np.max(self.data)

    def __add__(self, other):
        """Add two similar grids or a grid and a scaler: self + other."""
        if isinstance(other, SHGrid):
            if (self.grid == other.grid and self.data.shape ==
                    other.data.shape and self.kind == other.kind):
                data = self.data + other.data
                return SHGrid.from_array(data, grid=self.grid)
            else:
                raise ValueError('The two grids must be of the '
                                 'same kind and have the same shape and '
                                 'units.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not add a complex constant to a '
                                 'real grid.')
            data = self.data + other
            return SHGrid.from_array(data, grid=self.grid)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __radd__(self, other):
        """Add two similar grids or a grid and a scaler: self + other."""
        return self.__add__(other)

    def __sub__(self, other):
        """Subtract two similar grids or a grid and a scaler: self - other."""
        if isinstance(other, SHGrid):
            if (self.grid == other.grid and self.data.shape ==
                    other.data.shape and self.kind == other.kind):
                data = self.data - other.data
                return SHGrid.from_array(data, grid=self.grid)
            else:
                raise ValueError('The two grids must be of the '
                                 'same kind and have the same shape and '
                                 'units.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not subtract a complex constant from '
                                 'a real grid.')
            data = self.data - other
            return SHGrid.from_array(data, grid=self.grid)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __rsub__(self, other):
        """Subtract two similar grids or a grid and a scaler: other - self."""
        if isinstance(other, SHGrid):
            if (self.grid == other.grid and self.data.shape ==
                    other.data.shape and self.kind == other.kind):
                data = other.data - self.data
                return SHGrid.from_array(data, grid=self.grid)
            else:
                raise ValueError('The two grids must be of the '
                                 'same kind and have the same shape and '
                                 'units.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not subtract a complex constant from '
                                 'a real grid.')
            data = other - self.data
            return SHGrid.from_array(data, grid=self.grid)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __mul__(self, other):
        """Multiply two similar grids or a grid and a scaler: self * other."""
        if isinstance(other, SHGrid):
            if (self.grid == other.grid and self.data.shape ==
                    other.data.shape and self.kind == other.kind):
                data = self.data * other.data
                return SHGrid.from_array(data, grid=self.grid)
            else:
                raise ValueError('The two grids must be of the '
                                 'same kind and have the same shape.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not multiply a real grid by a complex '
                                 'constant.')
            data = self.data * other
            return SHGrid.from_array(data, grid=self.grid, units=self.units)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __rmul__(self, other):
        """Multiply two similar grids or a grid and a scaler: other * self."""
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        Divide two similar grids or a grid and a scalar.
        """
        if isinstance(other, SHGrid):
            if (self.grid == other.grid and self.data.shape ==
                    other.data.shape and self.kind == other.kind):
                data = self.data / other.data
                return SHGrid.from_array(data, grid=self.grid)
            else:
                raise ValueError('The two grids must be of the '
                                 'same kind and have the same shape.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not divide a real grid by a complex '
                                 'constant.')
            data = self.data / other
            return SHGrid.from_array(data, grid=self.grid, units=self.units)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __pow__(self, other):
        """Raise a grid to a scalar power: pow(self, other)."""
        if _np.isscalar(other) is True:
            return SHGrid.from_array(pow(self.data, other), grid=self.grid)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __abs__(self):
        """Return the absolute value of the gridded data."""
        return SHGrid.from_array(abs(self.data), grid=self.grid)

    def __repr__(self):
        str = ('kind = {:s}\n'
               'grid = {:s}\n'.format(repr(self.kind), repr(self.grid)))
        if self.grid == 'DH':
            str += ('n = {:d}\n'
                    'sampling = {:d}\n'.format(self.n, self.sampling))
        str += ('nlat = {:d}\n'
                'nlon = {:d}\n'
                'lmax = {:d}\n'
                'units = {:s}\n'
                'extend = {}'.format(self.nlat, self.nlon, self.lmax,
                                     repr(self.units), self.extend))
        return str

    # ---- Extract grid properties ----
    def lats(self, degrees=True):
        """
        Return the latitudes of each row of the gridded data.

        Usage
        -----
        lats = x.lats([degrees])

        Returns
        -------
        lats : ndarray, shape (nlat)
            1-D numpy array of size nlat containing the latitude of each row
            of the gridded data.

        Parameters
        ----------
        degrees : bool, optional, default = True
            If True, the output will be in degrees. If False, the output will
            be in radians.
        """
        if degrees is False:
            return _np.radians(self._lats())
        else:
            return self._lats()

    def lons(self, degrees=True):
        """
        Return the longitudes of each column of the gridded data.

        Usage
        -----
        lons = x.get_lon([degrees])

        Returns
        -------
        lons : ndarray, shape (nlon)
            1-D numpy array of size nlon containing the longitude of each row
            of the gridded data.

        Parameters
        ----------
        degrees : bool, optional, default = True
            If True, the output will be in degrees. If False, the output will
            be in radians.
        """
        if degrees is False:
            return _np.radians(self._lons())
        else:
            return self._lons()

    # ---- Functions that act on the data ----
    def histogram(self, bins=10, range=None):
        """
        Return an area-weighted histogram of the gridded data, normalized such
        that the integral over the range is unity.

        Usage
        -----
        hist, bin_edges = x.historgram([bins, range])

        Returns
        -------
        hist : array
            The values of the histogram, normalized such that the integral over
            the range in unity.
        bin_edges : array
            The values of the edges of the histogram bins.

        Parameters
        ----------
        bins : int or sequence of scalars or str, optional, default = 10
             If bins is an int, it defines the number of equal-width bins in
             the given range. If bins is a sequence, it defines a monotonically
             increasing array of bin edges, including the rightmost edge,
             allowing for non-uniform bin widths. If bins is a string, it
             defines the method used to calculate the optimal bin width, as
             defined by numpy.histogram_bin_edges.
        range : (float, float), optional, default = None
            The lower and upper range of the bins.

        Notes
        -----
        This method does not work with complex data.
        """
        return self._histogram(bins=bins, range=range)

    def expand(self, normalization='4pi', csphase=1, lmax_calc=None,
               backend=None, nthreads=None):
        """
        Expand the grid into spherical harmonics.

        Usage
        -----
        clm = x.expand([normalization, csphase, lmax_calc, backend, nthreads])

        Returns
        -------
        clm : SHCoeffs class instance

        Parameters
        ----------
        normalization : str, optional, default = '4pi'
            Normalization of the output class: '4pi', 'ortho', 'schmidt', or
            'unnorm', for geodesy 4pi normalized, orthonormalized, Schmidt
            semi-normalized, or unnormalized coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        lmax_calc : int, optional, default = x.lmax
            Maximum spherical harmonic degree to return.
        backend : str, optional, default = preferred_backend()
            Name of the preferred backend, either 'shtools' or 'ducc'.
        nthreads : int, optional, default = 1
            Number of threads to use for the 'ducc' backend. Setting this
            parameter to 0 will use as many threads as there are hardware
            threads on the system.

        Notes
        -----
        When expanding a Driscoll and Healy (1994) sampled grid (grid='DH' or
        'DH2') into spherical harmonic coefficients, the latitudinal bands at
        90 N and S are downweighted to zero and have no influence on the
        returned spherical harmonic coefficients.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        if backend is None:
            backend = preferred_backend()

        return self._expand(normalization=normalization, csphase=csphase,
                            lmax_calc=lmax_calc, backend=backend,
                            nthreads=nthreads)

    # ---- Plotting routines ----
    def plot3d(self, elevation=20, azimuth=30, cmap='viridis',
               cmap_limits=None, cmap_reverse=False, title=False,
               titlesize=None, scale=4., ax=None, show=True, fname=None):
        """
        Plot a 3-dimensional representation of the data.

        Usage
        -----
        x.plot3d([elevation, azimuth, cmap, cmap_limits, cmap_reverse, title,
                  titlesize, scale, ax, show, fname])

        Parameters
        ----------
        elevation : float, optional, default = 20
            The elevation angle (co-latitude), in degrees.
        azimuth : float, optional, default = 30
            The azimuth angle phi in degrees.
        cmap : str, optional, default = 'viridis'
            Name of the color map to use.
        cmap_reverse : bool, optional, default = False
            Set to True to reverse the sense of the color progression in the
            color table.
        cmap_limits : list, optional, default = [self.min(), self.max()]
            Set the lower and upper limits of the data used by the colormap,
            and optionally an interval for each color band. If the
            interval is specified, the number of discrete colors will be
            (cmap_limits[1]-cmap_limits[0])/cmap_limits[2].
        titlesize : int, optional, default = None
            The font size of the title.
        scale : float, optional, default = 4.
            The data (normalized by the maximum value) will be divided by scale
            before being added to the unit sphere.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear. This
            axes will be removed from the figure and replaced with an
            Axes3DSubplot object.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.

        Notes
        -----
        This 3-dimensional plotting routine plots the function

            1 + self.data / abs(self.data).max() / scale

        from a given elevation and azimuth viewing geometry. The data are first
        normalized by the maximum value, divived by scale, and then added to
        the unit sphere. If the data are complex, the magnitude of the data
        will be plotted. The color map corresponds to the values in self.data.

        This routine makes use of matplotlib3d and is very slow for large
        grids.
        """
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        if titlesize is None:
            titlesize = _mpl.rcParams['axes.titlesize']

        # make colormap
        if cmap_limits is None:
            cmap_limits = [self.min(), self.max()]
        if len(cmap_limits) == 3:
            num = int((cmap_limits[1] - cmap_limits[0]) / cmap_limits[2])
            if isinstance(cmap, _mpl.colors.Colormap):
                cmap_scaled = cmap._resample(num)
            else:
                cmap_scaled = _mpl.cm.get_cmap(cmap, num)
        else:
            cmap_scaled = _mpl.cm.get_cmap(cmap)
        if cmap_reverse:
            cmap_scaled = cmap_scaled.reversed()

        if ax is None:
            fig = _plt.figure()
            ax3d = fig.add_subplot(1, 1, 1, projection='3d')
        else:
            fig = ax.get_figure()
            subplotspec = ax.get_subplotspec()
            ax.remove()
            ax3d = fig.add_subplot(subplotspec, projection='3d')

        if self.kind == 'real':
            data = self.data
        elif self.kind == 'complex':
            data = _np.abs(self.data)
        else:
            raise ValueError('Grid has to be either real or complex, not {}.'
                             .format(self.kind))

        nlat, nlon = self.nlat, self.nlon
        lats = self.lats()
        lons = self.lons()

        if self.grid == 'DH':
            if self.extend:
                lats_circular = lats
            else:
                # add south pole
                lats_circular = _np.append(lats, [-90.])
        elif self.grid == 'GLQ':
            # add north and south pole
            lats_circular = _np.hstack(([90.], lats, [-90.]))
        if self.extend:
            lons_circular = lons
        else:
            lons_circular = _np.append(lons, [360.])

        nlats_circular = len(lats_circular)
        nlons_circular = len(lons_circular)
        sshape = nlats_circular, nlons_circular

        # make uv sphere and store all points
        phi = _np.radians(lons_circular)
        theta = _np.radians(90. - lats_circular)
        x = _np.sin(theta)[:, None] * _np.cos(phi)[None, :]
        y = _np.sin(theta)[:, None] * _np.sin(phi)[None, :]
        z = _np.cos(theta)[:, None] * _np.ones_like(lons_circular)[None, :]
        points = _np.vstack((x.flatten(), y.flatten(), z.flatten()))

        # The data need to be on a grid that spans 0-360 longitude and
        # -90-90 latitude. Add extra rows and columns if needed.
        if self.grid == 'DH':
            if self.extend:
                magn_point = data
            else:
                magn_point = _np.zeros((nlat + 1, nlon + 1))
                magn_point[:-1, :-1] = data
                magn_point[-1, :] = _np.mean(data[-1])
                magn_point[:-1, -1] = data[:, 0]
        if self.grid == 'GLQ':
            if self.extend:
                magn_point = _np.zeros((nlat + 2, nlon))
                magn_point[1:-1, :] = data
                magn_point[0, :] = _np.mean(data[0])
                magn_point[-1, :] = _np.mean(data[-1])
            else:
                magn_point = _np.zeros((nlat + 2, nlon + 1))
                magn_point[1:-1, :-1] = data
                magn_point[0, :] = _np.mean(data[0])
                magn_point[-1, :] = _np.mean(data[-1])
                magn_point[1:-1, -1] = data[:, 0]

        # compute face color, which is the average of all neighbour points
        magn_face = 1./4. * (magn_point[1:, 1:] + magn_point[:-1, 1:] +
                             magn_point[1:, :-1] + magn_point[:-1, :-1])

        magnmax_point = _np.max(_np.abs(magn_point))

        # compute colours and displace the points
        norm = _plt.Normalize(cmap_limits[0], cmap_limits[1])
        colors = cmap_scaled(norm(magn_face.flatten()))
        colors = colors.reshape(nlats_circular - 1, nlons_circular - 1, 4)
        points *= (1. + magn_point.flatten() / magnmax_point / scale)
        x = points[0].reshape(sshape)
        y = points[1].reshape(sshape)
        z = points[2].reshape(sshape)

        # plot data
        ax3d.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=colors,
                          shade=False)
        ax3d.set(xlim=(-1., 1.), ylim=(-1., 1.), zlim=(-1., 1.),
                 xticks=[-1, 1], yticks=[-1, 1], zticks=[-1, 1])
        ax3d.set_axis_off()
        ax3d.view_init(elev=elevation, azim=azimuth)

        if title:
            ax3d.set_title(title, fontsize=titlesize, pad=0.)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, ax3d

    def plot(self, projection=None, tick_interval=[30, 30], ticks='WSen',
             minor_tick_interval=[None, None], title=None, titlesize=None,
             colorbar=None, cmap='viridis', cmap_limits=None,
             cmap_limits_complex=None, cmap_reverse=False,
             cb_triangles='neither', cb_label=None, cb_ylabel=None,
             cb_tick_interval=None, cb_minor_tick_interval=None,
             cb_offset=None, cb_width=None, grid=False, axes_labelsize=None,
             tick_labelsize=None, xlabel=True, ylabel=True, ax=None, ax2=None,
             show=True, fname=None):
        """
        Plot the data using a Cartopy projection or a matplotlib cylindrical
        projection.

        Usage
        -----
        fig, ax = x.plot([projection, tick_interval, minor_tick_interval,
                          ticks, xlabel, ylabel, title, colorbar, cmap,
                          cmap_limits, cmap_limits_complex, cmap_reverse,
                          cb_triangles, cb_label, cb_ylabel, cb_tick_interval,
                          cb_minor_tick_interval, cb_offset, cb_width, grid,
                          titlesize, axes_labelsize, tick_labelsize, ax, ax2,
                          show, fname])

        Parameters
        ----------
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot, respectively. Alternatively, use
            'L', 'B', 'R', and 'T' for left, bottom, right, and top.
        xlabel : str, optional, default = 'Longitude' or 'GLQ longitude index'
            Label for the longitude axis.
        ylabel : str, optional, default = 'Latitude' or 'GLQ latitude index'
            Label for the latitude axis.
        title : str or list, optional, default = None
            The title of the plot. If the grid is complex, title should be a
            list of strings for the real and complex components.
        colorbar : str, optional, default = None
            Plot a colorbar along the 'top', 'right', 'bottom', or 'left' axis.
        cmap : str, optional, default = 'viridis'
            The color map to use when plotting the data and colorbar.
        cmap_limits : list, optional, default = [self.min(), self.max()]
            Set the lower and upper limits of the data used by the colormap,
            and optionally an interval for each color band. If the
            interval is specified, the number of discrete colors will be
            (cmap_limits[1]-cmap_limits[0])/cmap_limits[2]. If the data are
            complex, these limits will be used for the real component.
        cmap_limits_complex : list, optional, default = None
            Set the lower and upper limits of the imaginary component of the
            data used by the colormap, and optionally an interval for each
            color band.
        cmap_reverse : bool, optional, default = False
            Set to True to reverse the sense of the color progression in the
            color table.
        cb_triangles : str, optional, default = 'neither'
            Add triangles to the edges of the colorbar for minimum and maximum
            values. Can be 'neither', 'both', 'min', or 'max'.
        cb_label : str, optional, default = None
            Text label for the colorbar.
        cb_ylabel : str, optional, default = None
            Text label for the y axis of the colorbar
        cb_tick_interval : float, optional, default = None
            Colorbar major tick and annotation interval.
        cb_minor_tick_interval : float, optional, default = None
            Colorbar minor tick interval.
        cb_offset : float or int, optional, default = None
            Offset of the colorbar from the map edge in points. If None,
            the offset will be calculated automatically.
        cb_width : float, optional, default = None
            Width of the colorbar in percent with respect to the width of the
            respective image axis. Defaults are 2.5 and 5 for vertical and
            horizontal colorbars, respectively.
        grid : bool, optional, default = False
            If True, plot major grid lines.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        titlesize : int, optional, default = None
            The font size of the title.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear. If the
            grid is complex, the real component of the grid will be plotted
            on this axes.
        ax2 : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear. If the
            grid is complex, the complex component of the grid will be plotted
            on this axes.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
        """
        if projection is not None:
            if _cartopy_module:
                if not isinstance(projection, _ccrs.Projection):
                    raise ValueError('The input projection must be an '
                                     'instance of cartopy.crs.Projection.')
            else:
                raise ImportError('When using map projections, plot() '
                                  'requires installation of Cartopy.')
            if self.grid != 'DH':
                raise ValueError('Map projections are supported only for '
                                 'DH grid types.')

        if ticks is None:
            ticks = ''
        if tick_interval is None:
            tick_interval = [None, None]
        if minor_tick_interval is None:
            minor_tick_interval = [None, None]
        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
            if type(axes_labelsize) == str:
                axes_labelsize = _mpl.font_manager \
                                 .FontProperties(size=axes_labelsize) \
                                 .get_size_in_points()
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
            if type(tick_labelsize) == str:
                tick_labelsize = _mpl.font_manager \
                                 .FontProperties(size=tick_labelsize) \
                                 .get_size_in_points()
        if titlesize is None:
            titlesize = _mpl.rcParams['axes.titlesize']
            if type(titlesize) == str:
                titlesize = _mpl.font_manager \
                                 .FontProperties(size=titlesize) \
                                 .get_size_in_points()
        if self.kind == 'complex' and title is None:
            title = ['Real component', 'Imaginary component']
        if xlabel is True:
            if self.grid == 'DH':
                xlabel = 'Longitude'
            else:
                xlabel = 'GLQ longitude index'
        if ylabel is True:
            if self.grid == 'DH':
                ylabel = 'Latitude'
            else:
                ylabel = 'GLQ latitude index'
        if colorbar is not None:
            if colorbar not in set(['top', 'bottom', 'left', 'right']):
                raise ValueError("colorbar must be 'top', 'bottom', 'left' or "
                                 "'right'. Input value is {:s}."
                                 .format(repr(colorbar)))

        if ax is None and ax2 is None:
            fig, axes = self._plot(
                projection=projection, colorbar=colorbar,
                cb_triangles=cb_triangles, cb_label=cb_label, grid=grid,
                axes_labelsize=axes_labelsize, tick_labelsize=tick_labelsize,
                title=title, titlesize=titlesize, xlabel=xlabel, ylabel=ylabel,
                tick_interval=tick_interval, ticks=ticks,
                minor_tick_interval=minor_tick_interval,
                cb_tick_interval=cb_tick_interval, cb_ylabel=cb_ylabel,
                cb_minor_tick_interval=cb_minor_tick_interval, cmap=cmap,
                cmap_limits=cmap_limits, cb_offset=cb_offset,
                cb_width=cb_width, cmap_limits_complex=cmap_limits_complex,
                cmap_reverse=cmap_reverse)
        else:
            if self.kind == 'complex':
                if (ax is None and ax2 is not None) or (ax2 is None and
                                                        ax is not None):
                    raise ValueError('For complex grids, one must specify '
                                     'both optional arguments ax and ax2.')
            self._plot(projection=projection,
                       ax=ax, ax2=ax2, colorbar=colorbar,
                       cb_triangles=cb_triangles, cb_label=cb_label,
                       grid=grid, axes_labelsize=axes_labelsize,
                       tick_labelsize=tick_labelsize, title=title,
                       xlabel=xlabel, ylabel=ylabel,
                       tick_interval=tick_interval, ticks=ticks,
                       minor_tick_interval=minor_tick_interval,
                       titlesize=titlesize, cmap=cmap, cb_offset=cb_offset,
                       cb_tick_interval=cb_tick_interval,
                       cb_minor_tick_interval=cb_minor_tick_interval,
                       cmap_limits=cmap_limits, cb_ylabel=cb_ylabel,
                       cb_width=cb_width,
                       cmap_limits_complex=cmap_limits_complex,
                       cmap_reverse=cmap_reverse)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plotgmt(self, fig=None, projection='mollweide', region='g',
                width=None, unit='i', central_latitude=0, central_longitude=0,
                grid=[30, 30], tick_interval=[30, 30],
                minor_tick_interval=[None, None], ticks='WSen', title=None,
                cmap='viridis', cmap_limits=None, cmap_limits_complex=None,
                cmap_reverse=False, cmap_continuous=False, colorbar=None,
                cb_triangles='both', cb_label=None, cb_ylabel=None,
                cb_tick_interval=None, cb_minor_tick_interval=None,
                cb_offset=None, shading=None, shading_azimuth=-45.,
                shading_amplitude=1.0, titlesize=None, axes_labelsize=None,
                tick_labelsize=None, horizon=60, offset=[None, None],
                fname=None):
        """
        Plot projected data using the Generic Mapping Tools (pygmt).

        To display the figure in a jupyter notebook, use
            fig.show()
        To display the figure in the terminal environment, use
            fig.show(method='external')

        Usage
        -----
        fig = x.plotgmt([fig, projection, region, width, unit,
                         central_latitude, central_longitude, grid,
                         tick_interval, minor_tick_interval, ticks, title,
                         cmap, cmap_limits, cmap_limits_complex, cmap_reverse,
                         cmap_continuous, colorbar, cb_triangles, cb_label,
                         cb_ylabel, cb_tick_interval, cb_minor_tick_interval,
                         cb_offset, shading, shading_azimuth,
                         shading_amplitude, titlesize, axes_labelsize,
                         tick_labelsize, horizon, offset, fname])

        Returns
        -------
        fig : pygmt.figure.Figure class instance

        Parameters
        ----------
        fig : pygmt.Figure() class instance, optional, default = None
            If provided, the plot will be placed in a pre-existing figure.
        projection : str, optional, default = 'mollweide'
            The name of a global or hemispherical projection (see Notes). Only
            the first three characters are necessary to identify the
            projection.
        region : str or list, optional, default = 'g'
            The map region, consisting of a list [West, East, South, North] in
            degrees. The default 'g' specifies the entire sphere.
        width : float, optional, default = mpl.rcParams['figure.figsize'][0]
            The width of the projected image.
        unit : str, optional, default = 'i'
            The measurement unit of the figure width and offset: 'i' for
            inches or 'c' for cm.
        central_longitude : float, optional, default = 0
            The central meridian or center of the projection.
        central_latitude : float, optional, default = 0
            The center of the projection used with hemispheric projections, or
            the standard parallel used with cylindrical projections.
        grid : list, optional, default = [30, 30]
            Grid line interval [longitude, latitude] in degrees. If None, grid
            lines will not be plotted for that axis. If true, gridlines will be
            plotted at the same interval as the major ticks.
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters will plot the ticks and annotations, whereas small letters
            will plot only the ticks. 'W', 'S', 'E', and 'N' denote the west,
            south, east and north boundaries of the plot.
        title : str, optional, default = None
            The title to be displayed above the plot.
        cmap : str, optional, default = 'viridis'
            The color map to use when plotting the data and colorbar.
        cmap_limits : list, optional, default = [self.min(), self.max()]
            Set the lower and upper limits of the data used by the colormap
            when plotting, and optionally an interval for each color.
        cmap_limits_complex : list, optional, default = None
            Set the lower and upper limits of the imaginary component of the
            data used by the colormap, and optionally an interval for each
            color band.
        cmap_reverse : bool, optional, default = False
            Set to True to reverse the sense of the color progression in the
            color table.
        cmap_continuous : bool, optional, default = False
            If True, create a continuous colormap. Default behavior is to
            use contant colors for each interval.
        colorbar : str, optional, default = None
            Plot a colorbar along the 'top', 'right', 'bottom', or 'left' axis.
        cb_triangles : str, optional, default = 'both'
            Add triangles to the edges of the colorbar for minimum and maximum
            values. Can be 'neither', 'both', 'min', or 'max'.
        cb_label : str, optional, default = None
            Text label for the colorbar.
        cb_ylabel : str, optional, default = None
            Text label for the y axis of the colorbar
        cb_tick_interval : float, optional, default = None
            Annotation interval on the colorbar.
        cb_minor_tick_interval : float, optional, default = None
            Colorbar minor tick interval.
        cb_offset : float or int, optional, default = None
            Offset of the colorbar from the map edge in points. If None,
            the offset will be calculated automatically.
        shading : bool, str, or SHGrid instance, optional, default = None
            Apply intensity shading to the image. The shading (with values
            from -1 to 1) can be derived from the data by setting to True,
            from an external netcdf file by supplying a filename, or from an
            SHGrid class instance by supplying the name of the SHGrid. When
            intensity shading is applied, the default behavior is to create a
            gradient of the shading data. If it is not necessary to create a
            gradient, shading_azimuth should be set to None. If shading is
            None, no intensity shading will be applied.
        shading_azimuth : float, optional, default = -45.
            When applying intensity shading to the image, a gradient of the
            shading data is computed using the supplied azimuth direction (in
            degrees). If it is not necessary to create a gradient from the
            shading data, shading_azimuth should be set to None. When shading
            is set to True, a shading azimuth must be provided.
        shading_amplitude : float, optional, default = 1.
            The maximum amplitude of the intensity used in the shading, from
            0 to 1.
        titlesize : int, optional, default = None
            The font size of the title.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        horizon : float, optional, default = 60
            The horizon (number of degrees from the center to the edge) used
            with the Gnomonic projection.
        offset : list, optional, default = [None, None]
            Offset of the plot in the x and y directions from the current
            origin.
        fname : str, optional, default = None
            If present, save the image to the specified file.

        Notes
        -----
        Global and hemispherical projections (region='g') with corresponding
        abbreviation used by `projection`:

        Azimuthal projections
            Lambert-azimuthal-equal-area (lam)
            Stereographic-equal-angle (ste)
            Orthographic (ort)
            Azimuthal-equidistant (azi)
            Gnomonic (gno)

        Cylindrical projections (case sensitive)
            cylindrical-equidistant (cyl)
            Cylindrical-equal-area (Cyl)
            CYLindrical-stereographic (CYL)
            Miller-cylindrical (mil)

        Miscellaneous projections
            Mollweide (mol)
            Hammer (ham)
            Winkel-Tripel (win)
            Robinson (rob)
            Eckert (eck)
            Sinusoidal (sin)
            Van-der-Grinten (van)
        """
        if not _pygmt_module:
            raise ImportError('plotgmt() requires installation of the module '
                              'pygmt.')
        if tick_interval is None:
            tick_interval = [None, None]
        if minor_tick_interval is None:
            minor_tick_interval = [None, None]
        if self.kind == 'complex' and title is None:
            title = ['Real component', 'Imaginary component']
        if grid is True:
            grid = tick_interval
        if width is None:
            width = _mpl.rcParams['figure.figsize'][0]
        if colorbar is not None:
            if colorbar not in set(['top', 'bottom', 'left', 'right']):
                raise ValueError("colorbar must be 'top', 'bottom', 'left' or "
                                 "'right'. Input value is {:s}."
                                 .format(repr(colorbar)))
        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
            if type(axes_labelsize) == str:
                axes_labelsize = _mpl.font_manager \
                                 .FontProperties(size=axes_labelsize) \
                                 .get_size_in_points()
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
            if type(tick_labelsize) == str:
                tick_labelsize = _mpl.font_manager \
                                 .FontProperties(size=tick_labelsize) \
                                 .get_size_in_points()
        if titlesize is None:
            titlesize = _mpl.rcParams['axes.titlesize']
            if type(titlesize) == str:
                titlesize = _mpl.font_manager \
                                 .FontProperties(size=titlesize) \
                                 .get_size_in_points()

        figure = self._plot_pygmt(
            fig=fig, projection=projection, region=region, width=width,
            unit=unit, central_latitude=central_latitude,
            central_longitude=central_longitude, grid=grid,
            tick_interval=tick_interval,
            minor_tick_interval=minor_tick_interval, ticks=ticks, title=title,
            cmap=cmap, cmap_limits=cmap_limits,
            cmap_limits_complex=cmap_limits_complex, cmap_reverse=cmap_reverse,
            cmap_continuous=cmap_continuous, colorbar=colorbar,
            cb_triangles=cb_triangles, cb_label=cb_label, cb_ylabel=cb_ylabel,
            cb_tick_interval=cb_tick_interval, cb_offset=cb_offset,
            cb_minor_tick_interval=cb_minor_tick_interval, shading=shading,
            shading_azimuth=shading_azimuth,
            shading_amplitude=shading_amplitude, titlesize=titlesize,
            axes_labelsize=axes_labelsize, tick_labelsize=tick_labelsize,
            horizon=horizon, offset=offset)

        if fname is not None:
            figure.savefig(fname)
        if fig is None:
            return figure

    def plot_histogram(self, bins=10, range=None, cumulative=False,
                       histtype='bar', orientation='vertical', xscale='lin',
                       yscale='lin', color=None, legend=None,
                       legend_loc='best', xlabel=None,
                       ylabel='Fraction', title=None, grid=False,
                       axes_labelsize=None, tick_labelsize=None,
                       titlesize=None, ax=None, show=True, fname=None):
        """
        Plot an area-weighted histogram of the gridded data, normalized such
        that the integral over the range is unity.

        Usage
        -----
        x.plot_histogram([bins, range, cumulative, histtype, orientation,
                          xscale, yscale, color, legend, legend_loc, xlabel,
                          ylabel, title, grid, axes_labelsize, tick_labelsize,
                          titlesize, ax, show, fname])

        Parameters
        ----------
        bins : int or sequence of scalars or str, optional, default = 10
             If bins is an int, it defines the number of equal-width bins in
             the given range. If bins is a sequence, it defines a monotonically
             increasing array of bin edges, including the rightmost edge,
             allowing for non-uniform bin widths. If bins is a string, it
             defines the method used to calculate the optimal bin width, as
             defined by numpy.histogram_bin_edges.
        range : (float, float), optional, default = None
            The lower and upper range of the bins.
        cumulative : bool or -1, optional, default = False
            If True, then a histogram is computed where each bin gives the
            counts in that bin plus all bins for smaller values. If cumulative
            is -1, the direction of accumulation is reversed.
        histtype : str, optional, default = 'bar'
            The type of histogram to draw. 'bar' is a traditional bar-type
            histogram. 'barstacked' is a bar-type histogram where multiple data
            are stacked on top of each other. 'step' generates a lineplot that
            is unfilled. 'stepfilled' generates a lineplot that is filled.
        orientation : str, optional, default = 'vertical'
            Orientiation of the histogram, either 'vertical' or 'horizontal'.
        xscale : str, optional, default = 'lin'
            Scale of the x axis: 'lin' for linear or 'log' for logarithmic.
        yscale : str, optional, default = 'lin'
            Scale of the y axis: 'lin' for linear or 'log' for logarithmic.
        title : str or list, optional, default = None
            The title of the plot.
        color : str, optional, default = None
            Name of the color used for the histogram bars or lines.
        legend : str or None, optional, default = None
            Label to use when plotting multiple datasets.
        legend_loc : str, optional, default = 'best'
            Location of the legend, such as 'upper right' or 'lower center'
            (see pyplot.legend for all options).
        xlabel : str, optional, default = None
            Label for the x axis.
        ylabel : str, optional, default = 'Fraction'
            Label for the y axis.
        grid : bool, optional, default = False
            If True, plot grid lines.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        titlesize : int, optional, default = None
            The font size of the title.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.

        Notes
        -----
        This method calls histogram() to bin the data, and then uses
        matplotlib.pyplot.hist() to plot the binned data. This method does not
        work with complex data.
        """
        if self.kind == 'complex':
            raise NotImplementedError(
                'plot_histogram() does not support complex data.')

        if ax is None:
            fig, axes = _plt.subplots(1, 1)
        else:
            axes = ax

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
            if type(axes_labelsize) == str:
                axes_labelsize = _mpl.font_manager \
                                 .FontProperties(size=axes_labelsize) \
                                 .get_size_in_points()
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
            if type(tick_labelsize) == str:
                tick_labelsize = _mpl.font_manager \
                                 .FontProperties(size=tick_labelsize) \
                                 .get_size_in_points()
        if titlesize is None:
            titlesize = _mpl.rcParams['axes.titlesize']
            if type(titlesize) == str:
                titlesize = _mpl.font_manager \
                                 .FontProperties(size=titlesize) \
                                 .get_size_in_points()

        hist, bin_edges = self.histogram(bins=bins, range=range)
        axes.hist(bin_edges[:-1], bin_edges, weights=hist,
                  cumulative=cumulative, histtype=histtype,
                  orientation=orientation, color=color, label=legend)

        if xscale == 'log':
            axes.set_xscale('log')
        if yscale == 'log':
            axes.set_yscale('log')

        if xlabel is not None:
            axes.set_xlabel(xlabel, fontsize=axes_labelsize)
        if ylabel is not None:
            axes.set_ylabel(ylabel, fontsize=axes_labelsize)
        if title is not None:
            axes.set_title(title, fontsize=titlesize)

        axes.tick_params(which='major', labelsize=tick_labelsize)
        axes.minorticks_on()
        axes.grid(grid, which='major')

        if legend is not None:
            axes.legend(loc=legend_loc)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes


# ---- Real Driscoll and Healy grid class ----

class DHRealGrid(SHGrid):
    """Class for real Driscoll and Healy (1994) grids."""

    @staticmethod
    def istype(kind):
        return kind == 'real'

    @staticmethod
    def isgrid(grid):
        return grid == 'DH'

    def __init__(self, array, units=None, copy=True):
        self.nlat, self.nlon = array.shape

        if self.nlat % 2 != 0:
            self.n = self.nlat - 1
            self.extend = True
        else:
            self.n = self.nlat
            self.extend = False

        if self.nlon == 2 * self.nlat - self.extend:
            self.sampling = 2
        elif self.nlon == self.nlat:
            self.sampling = 1
        else:
            raise ValueError('Input array has shape (nlat={:d}, nlon={:d}) '
                             .format(self.nlat, self.nlon) +
                             'but needs nlon=nlat, nlon=2*nlat, or '
                             'nlon=2*nlat-1.'
                             )

        self.lmax = int(self.n / 2 - 1)
        self.grid = 'DH'
        self.kind = 'real'
        self.units = units

        if copy:
            self.data = _np.copy(array)
        else:
            self.data = array

    def _lats(self):
        """Return the latitudes (in degrees) of the gridded data."""
        if self.extend:
            lats = _np.linspace(90.0, -90.0, num=self.nlat)
        else:
            lats = _np.linspace(90.0, -90.0 + 180.0 / self.nlat, num=self.nlat)
        return lats

    def _lons(self):
        """Return the longitudes (in degrees) of the gridded data."""
        if self.extend:
            lons = _np.linspace(0.0, 360.0, num=self.nlon)
        else:
            lons = _np.linspace(0.0, 360.0 - 360.0 / self.nlon, num=self.nlon)
        return lons

    def _histogram(self, bins=None, range=None):
        """Return an area-weighted histogram normalized to unity."""
        delta_phi = self.lons()[1] - self.lons()[0]
        delta_phi *= _np.pi / 180.

        da = _np.zeros_like(self.data)

        # i=0, 90 N
        theta2 = 90.0 - (90.0 + self.lats()[1]) / 2.
        da[0, :] = 1. - _np.cos(theta2 * _np.pi / 180.)
        for i in _np.arange(1, self.nlat-1):
            theta1 = 90.0 - (self.lats()[i-1] + self.lats()[i]) / 2.
            theta2 = 90.0 - (self.lats()[i] + self.lats()[i+1]) / 2.
            da[i, :] = _np.cos(theta1 * _np.pi / 180.) - \
                _np.cos(theta2 * _np.pi / 180.)
        # last latitudinal band
        i = self.nlat - 1
        if self.extend:
            theta1 = 90.0 - (-90.0 + self.lats()[i-1]) / 2.
            da[i, :] = _np.cos(theta1 * _np.pi / 180.) + 1
        else:
            theta1 = 90.0 - (self.lats()[i-1] + self.lats()[i]) / 2.
            theta2 = 90.0 - (self.lats()[i] - 90.) / 2.
            da[i, :] = _np.cos(theta1 * _np.pi / 180.) - \
                _np.cos(theta2 * _np.pi / 180.)

        return _np.histogram(self.data[:, :self.nlon-self.extend],
                             bins=bins,
                             weights=da[:, :self.nlon-self.extend],
                             density=True, range=range)

    def _expand(self, normalization, csphase, lmax_calc, backend, nthreads):
        """Expand the grid into real spherical harmonics."""
        from .shcoeffs import SHCoeffs
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'unnorm':
            norm = 3
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        cilm = backend_module(
            backend=backend, nthreads=nthreads).SHExpandDH(
                self.data[:self.nlat-self.extend, :self.nlon-self.extend],
                norm=norm, csphase=csphase, sampling=self.sampling,
                lmax_calc=lmax_calc)
        coeffs = SHCoeffs.from_array(cilm,
                                     normalization=normalization.lower(),
                                     csphase=csphase, units=self.units,
                                     copy=False)
        return coeffs

    def _plot(self, projection=None, xlabel=None, ylabel=None, ax=None,
              ax2=None, colorbar=None, cb_triangles=None, cb_label=None,
              grid=False, axes_labelsize=None, tick_labelsize=None, title=None,
              titlesize=None, cmap=None, tick_interval=None, ticks=None,
              minor_tick_interval=None, cb_tick_interval=None, cb_ylabel=None,
              cb_minor_tick_interval=None, cmap_limits=None, cmap_reverse=None,
              cmap_limits_complex=None, cb_offset=None, cb_width=None):
        """Plot the data as a matplotlib cylindrical projection,
           or with Cartopy when projection is specified."""
        if ax is None:
            if colorbar is not None:
                if colorbar in set(['top', 'bottom']):
                    scale = 0.67
                else:
                    scale = 0.5
            else:
                scale = 0.55
            figsize = (_mpl.rcParams['figure.figsize'][0],
                       _mpl.rcParams['figure.figsize'][0] * scale)
            fig = _plt.figure(figsize=figsize)
            if projection is not None:
                axes = fig.add_subplot(111, projection=projection)
            else:
                axes = fig.add_subplot(111)
        else:
            if projection is not None:
                fig = ax.get_figure()
                subplotspec = ax.get_subplotspec()
                ax.remove()
                axes = fig.add_subplot(subplotspec, projection=projection)
            else:
                axes = ax

        # set tick intervals
        if tick_interval[0] is None:
            xticks = []
        else:
            xticks = _np.linspace(0, 360, num=360//tick_interval[0]+1,
                                  endpoint=True)
        if tick_interval[1] is None:
            yticks = []
        else:
            yticks = _np.linspace(-90, 90, num=180//tick_interval[1]+1,
                                  endpoint=True)
        if minor_tick_interval[0] is None:
            minor_xticks = []
        else:
            minor_xticks = _np.linspace(
                0, 360, num=360//minor_tick_interval[0]+1, endpoint=True)
        if minor_tick_interval[1] is None:
            minor_yticks = []
        else:
            minor_yticks = _np.linspace(
                -90, 90, num=180//minor_tick_interval[1]+1, endpoint=True)

        # make colormap
        if cmap_limits is None:
            cmap_limits = [self.min(), self.max()]
        if len(cmap_limits) == 3:
            num = int((cmap_limits[1] - cmap_limits[0]) / cmap_limits[2])
            if isinstance(cmap, _mpl.colors.Colormap):
                cmap_scaled = cmap._resample(num)
            else:
                cmap_scaled = _mpl.cm.get_cmap(cmap, num)
        else:
            cmap_scaled = _mpl.cm.get_cmap(cmap)
        if cmap_reverse:
            cmap_scaled = cmap_scaled.reversed()

        # compute colorbar ticks
        cb_ticks = None
        cb_minor_ticks = None
        vmin = cmap_limits[0]
        vmax = cmap_limits[1]
        if cb_tick_interval is not None:
            if _np.sign(vmin) == -1.:
                start = (abs(vmin) // cb_tick_interval) \
                    * cb_tick_interval * _np.sign(vmin)
            else:
                start = (vmin // cb_tick_interval + 1) \
                    * cb_tick_interval
            if _np.sign(vmax) == -1.:
                stop = (abs(vmax) // cb_tick_interval + 1) \
                    * cb_tick_interval * _np.sign(vmax)
            else:
                stop = (vmax // cb_tick_interval) * cb_tick_interval
            cb_ticks = _np.linspace(start, stop,
                                    num=int((stop-start)//cb_tick_interval+1),
                                    endpoint=True)
        if cb_minor_tick_interval is not None:
            if _np.sign(vmin) == -1.:
                start = (abs(vmin) // cb_minor_tick_interval) \
                    * cb_minor_tick_interval * _np.sign(vmin)
            else:
                start = (vmin // cb_minor_tick_interval + 1) \
                    * cb_minor_tick_interval
            if _np.sign(vmax) == -1.:
                stop = (abs(vmax) // cb_minor_tick_interval + 1) \
                    * cb_minor_tick_interval * _np.sign(vmax)
            else:
                stop = (vmax // cb_minor_tick_interval) * \
                    cb_minor_tick_interval
            cb_minor_ticks = _np.linspace(
                start, stop, num=int((stop-start)//cb_minor_tick_interval+1),
                endpoint=True)

        # determine which ticks to plot
        if 'W' in ticks or 'L' in ticks:
            left, labelleft = True, True
        elif 'w' in ticks or 'l' in ticks:
            left, labelleft = True, False
        else:
            left, labelleft = False, False
        if 'S' in ticks or 'B' in ticks:
            bottom, labelbottom = True, True
        elif 's' in ticks or 'b' in ticks:
            bottom, labelbottom = True, False
        else:
            bottom, labelbottom = False, False
        if 'E' in ticks or 'R' in ticks:
            right, labelright = True, True
        elif 'e' in ticks or 'r' in ticks:
            right, labelright = True, False
        else:
            right, labelright = False, False
        if 'N' in ticks or 'T' in ticks:
            top, labeltop = True, True
        elif 'n' in ticks or 't' in ticks:
            top, labeltop = True, False
        else:
            top, labeltop = False, False

        extent = (-360. / self.sampling / self.n / 2.,
                  360. + 360. / self.sampling / self.n * (self.extend - 0.5),
                  -90. - 180. / self.n * (self.extend - 0.5),
                  90. + 180. / 2. / self.n)

        # Add space for annotations between plot and colorbar. This will be
        # False for map projections that do not support longitude labels
        cb_space = True

        # plot image, ticks, and annotations
        if projection is not None:
            axes.set_global()
            cim = axes.imshow(
                self.data, transform=_ccrs.PlateCarree(central_longitude=0.0),
                origin='upper', extent=extent, cmap=cmap_scaled,
                vmin=cmap_limits[0], vmax=cmap_limits[1])
            if isinstance(projection, _ccrs.PlateCarree):
                axes.set_xticks(
                    xticks, crs=_ccrs.PlateCarree(central_longitude=0.0))
                axes.set_yticks(
                    yticks, crs=_ccrs.PlateCarree(central_longitude=0.0))
                axes.set_xticks(minor_xticks, minor=True,
                                crs=_ccrs.PlateCarree(central_longitude=0.0))
                axes.set_yticks(minor_yticks, minor=True,
                                crs=_ccrs.PlateCarree(central_longitude=0.0))
                axes.xaxis.set_major_formatter(_LongitudeFormatter())
                axes.yaxis.set_major_formatter(_LatitudeFormatter())
            else:
                cb_space = False
            if grid:
                axes.gridlines(xlocs=xticks-180, ylocs=yticks,
                               crs=_ccrs.PlateCarree(central_longitude=0.0))
        else:
            cim = axes.imshow(self.data, origin='upper', extent=extent,
                              cmap=cmap_scaled, vmin=cmap_limits[0],
                              vmax=cmap_limits[1])
            axes.set(xlim=(0, 360), ylim=(-90, 90))
            axes.set_xlabel(xlabel, fontsize=axes_labelsize)
            axes.set_ylabel(ylabel, fontsize=axes_labelsize)
            axes.set(xticks=xticks, yticks=yticks)
            axes.set_xticks(minor_xticks, minor=True)
            axes.set_yticks(minor_yticks, minor=True)
            axes.grid(grid, which='major')
            axes.xaxis.set_major_formatter(
                _mpl.ticker.FormatStrFormatter("%d"+u"\u00B0"))
            axes.yaxis.set_major_formatter(
                _mpl.ticker.FormatStrFormatter("%d"+u"\u00B0"))

        axes.tick_params(bottom=bottom, top=top, right=right, left=left,
                         labelbottom=labelbottom, labeltop=labeltop,
                         labelleft=labelleft, labelright=labelright,
                         which='both')
        axes.tick_params(which='major', labelsize=tick_labelsize)
        if title is not None:
            axes.set_title(title, fontsize=titlesize)

        # plot colorbar
        if colorbar is not None:
            if cb_offset is None:
                if colorbar in set(['left', 'right']):
                    offset = 0.15
                    if (colorbar == 'left' and
                        ('W' in ticks or 'L' in ticks)) or \
                            (colorbar == 'right' and
                             ('E' in ticks or 'R' in ticks)):
                        offset += 2 * tick_labelsize / 72.
                    # add space for ylabel on left of plot only
                    if ylabel != '' and ylabel is not None and \
                            projection is None and colorbar == 'left':
                        offset += 1.4 * axes_labelsize / 72.
                else:
                    offset = 0.
                    # add space for ticks
                    if (colorbar == 'bottom' and bottom and cb_space) or \
                            (colorbar == 'top' and top and cb_space):
                        offset += _mpl.rcParams['xtick.major.size']
                    # add space for labels
                    if (colorbar == 'bottom' and labelbottom and cb_space) or \
                            (colorbar == 'top' and labeltop and cb_space):
                        offset += _mpl.rcParams['xtick.major.pad']
                        offset += tick_labelsize
                    # add space for xlabel on bottom of plot only
                    if xlabel != '' and xlabel is not None and \
                            projection is None and colorbar == 'bottom':
                        offset += axes_labelsize
                    offset += 1.3 * _mpl.rcParams['font.size']  # add extra
                    offset /= 72.  # convert to inches
            else:
                offset = cb_offset / 72.0  # convert to inches

            divider = _make_axes_locatable(axes)
            if colorbar in set(['left', 'right']):
                orientation = 'vertical'
                extendfrac = 0.05
                if cb_width is None:
                    size = '2.5%'
                else:
                    size = '{:f}%'.format(cb_width)
            else:
                orientation = 'horizontal'
                extendfrac = 0.025
                if cb_width is None:
                    size = '5%'
                else:
                    size = '{:f}%'.format(cb_width)
            cax = divider.append_axes(colorbar, size=size, pad=offset,
                                      axes_class=_plt.Axes)
            cbar = _plt.colorbar(cim, cax=cax, orientation=orientation,
                                 extend=cb_triangles, extendfrac=extendfrac)
            if colorbar == 'left':
                cbar.ax.yaxis.set_ticks_position('left')
                cbar.ax.yaxis.set_label_position('left')
            if colorbar == 'top':
                cbar.ax.xaxis.set_ticks_position('top')
                cbar.ax.xaxis.set_label_position('top')
            if cb_label is not None:
                cbar.set_label(cb_label, fontsize=axes_labelsize)
            if cb_ylabel is not None:
                if colorbar in set(['left', 'right']):
                    cbar.ax.xaxis.set_label_position('top')
                    cbar.ax.set_xlabel(cb_ylabel, fontsize=tick_labelsize)
                else:
                    cbar.ax.yaxis.set_label_position('right')
                    cbar.ax.set_ylabel(cb_ylabel, fontsize=tick_labelsize,
                                       rotation=0., labelpad=axes_labelsize/2.,
                                       va='center', ha='left')
            if cb_ticks is not None:
                cbar.set_ticks(cb_ticks)
            if cb_minor_ticks is not None:
                if colorbar in set(['top', 'bottom']):
                    cbar.ax.xaxis.set_ticks(cb_minor_ticks, minor=True)
                else:
                    cbar.ax.yaxis.set_ticks(cb_minor_ticks, minor=True)
            cbar.ax.tick_params(labelsize=tick_labelsize)

        if ax is None:
            return fig, axes

    def _plot_pygmt(self, fig=None, projection=None, region=None, width=None,
                    unit=None, central_latitude=None, central_longitude=None,
                    grid=None, tick_interval=None, minor_tick_interval=None,
                    ticks=None, title=None, cmap=None, cmap_limits=None,
                    cmap_limits_complex=None, cmap_reverse=None,
                    cmap_continuous=None, colorbar=None, cb_triangles=None,
                    cb_label=None, cb_ylabel=None, cb_tick_interval=None,
                    cb_minor_tick_interval=None, shading=None,
                    shading_azimuth=None, shading_amplitude=None,
                    titlesize=None, axes_labelsize=None, tick_labelsize=None,
                    horizon=None, offset=[None, None], cb_offset=None):
        """
        Plot projected data using pygmt.
        """
        center = [central_longitude, central_latitude]

        if projection.lower()[0:3] == 'mollweide'[0:3]:
            proj_str = 'W' + str(center[0])
        elif projection.lower()[0:3] == 'hammer'[0:3]:
            proj_str = 'H' + str(center[0])
        elif projection.lower()[0:3] == 'winkel-tripel'[0:3]:
            proj_str = 'R' + str(center[0])
        elif projection.lower()[0:3] == 'robinson'[0:3]:
            proj_str = 'N' + str(center[0])
        elif projection.lower()[0:3] == 'eckert'[0:3]:
            proj_str = 'K' + str(center[0])
        elif projection.lower()[0:3] == 'sinusoidal'[0:3]:
            proj_str = 'I' + str(center[0])
        elif projection.lower()[0:3] == 'van-der-grinten'[0:3]:
            proj_str = 'V' + str(center[0])
        elif projection.lower()[0:3] == 'lambert-azimuthal-equal-area'[0:3]:
            proj_str = 'A' + str(center[0]) + '/' + str(center[1])
        elif projection.lower()[0:3] == 'stereographic-equal-angle'[0:3]:
            proj_str = 'S' + str(center[0]) + '/' + str(center[1])
        elif projection.lower()[0:3] == 'orthographic'[0:3]:
            proj_str = 'G' + str(center[0]) + '/' + str(center[1])
        elif projection.lower()[0:3] == 'azimuthal-equidistant'[0:3]:
            proj_str = 'E' + str(center[0]) + '/' + str(center[1])
        elif projection.lower()[0:3] == 'gnomonic'[0:3]:
            proj_str = 'F' + str(center[0]) + '/' + str(center[1]) + '/' \
                + str(horizon)
        elif projection.lower()[0:3] == 'miller-cylindrical'[0:3]:
            proj_str = 'J' + str(central_longitude)
        elif projection[0:3] == 'cylindrical-equidistant'[0:3]:
            proj_str = 'Q' + str(center[0]) + '/' + str(center[1])
        elif projection[0:3] == 'Cylindrical-equal-area'[0:3]:
            proj_str = 'Y' + str(center[0]) + '/' + str(center[1])
        elif projection[0:3] == 'CYLindrical-stereographic'[0:3]:
            proj_str = 'Cyl_stere' + '/' + str(center[0]) + '/' \
                + str(center[1])
        else:
            raise ValueError('Input projection is not recognized or '
                             'supported. Input projection = {:s}'
                             .format(repr(projection)))

        proj_str += '/' + str(width) + unit

        framex = 'x'
        framey = 'y'
        if grid[0] is not None:
            framex += 'g' + str(grid[0])
        if grid[1] is not None:
            framey += 'g' + str(grid[1])
        if tick_interval[0] is not None:
            framex += 'a' + str(tick_interval[0])
        if tick_interval[1] is not None:
            framey += 'a' + str(tick_interval[1])
        if minor_tick_interval[0] is not None:
            framex += 'f' + str(minor_tick_interval[0])
        if minor_tick_interval[1] is not None:
            framey += 'f' + str(minor_tick_interval[1])
        if title is not None:
            ticks += '+t"{:s}"'.format(title)
        frame = [framex, framey, ticks]

        position = None
        cb_str = None
        if colorbar is not None:
            if colorbar == 'right':
                position = "JMR"
            elif colorbar == 'bottom':
                position = "JBC+h"
            elif colorbar == 'left':
                position = "JML"
            elif colorbar == 'top':
                position = "JTC+h"
            if cb_offset is not None:
                if colorbar == 'bottom' or colorbar == 'top':
                    position += '+o0p/' + str(cb_offset) + 'p'
                else:
                    position += '+o' + str(cb_offset) + 'p/0p'
            if cb_triangles is not None:
                if cb_triangles == 'neither':
                    pass
                elif cb_triangles == 'both':
                    position += '+ebf'
                elif cb_triangles == 'min':
                    position += '+eb'
                elif cb_triangles == 'max':
                    position += '+ef'
                else:
                    raise ValueError("cb_triangles must be 'neither', 'both' "
                                     "'min' or 'max'. Input value is {:s}."
                                     .format(repr(cb_triangles)))
            cb_str = []
            x_str = 'x'
            if cb_tick_interval is not None:
                x_str += 'a' + str(cb_tick_interval)
            if cb_minor_tick_interval is not None:
                x_str += 'f' + str(cb_minor_tick_interval)
            cb_str.extend([x_str])
            if cb_label is not None:
                cb_str.extend(['x+l"{:s}"'.format(cb_label)])
            if cb_ylabel is not None:
                cb_str.extend(['y+l"{:s}"'.format(cb_ylabel)])

        if offset[0] is not None:
            xshift = str(offset[0]) + unit
        else:
            xshift = None
        if offset[1] is not None:
            yshift = str(offset[1]) + unit
        else:
            yshift = None

        if fig is None:
            figure = _pygmt.Figure()
        else:
            figure = fig

        if cmap_limits is None:
            cmap_limits = [self.min(), self.max()]

        if shading is True:
            shading_str = "+a{:}+nt{:}+m0".format(shading_azimuth,
                                                  shading_amplitude)
        elif type(shading) is str:
            shading_str = shading
            if shading_azimuth is not None:
                shading_str += "+a{:}+nt{:}+m0".format(shading_azimuth,
                                                       shading_amplitude)
        elif isinstance(shading, SHGrid):
            if self.data.shape != shading.data.shape:
                raise ValueError('The input SHGrid used for shading '
                                 'must have the same shape as the grid being '
                                 'plotted. Shape of grid = {:}. '
                                 .format(self.data.shape) +
                                 'Shape of shading grid = {:}.'
                                 .format(shading.data.shape)
                                 )

            f = _tempfile.NamedTemporaryFile(prefix='shtools_', suffix='.nc')
            shading.to_netcdf(f.name)
            shading_str = f.name
            if shading_azimuth is not None:
                shading_str += "+a{:}+nt{:}+m0".format(shading_azimuth,
                                                       shading_amplitude)
        else:
            shading_str = None

        # Necessary to fix bug in pygmt 0.4+
        if shading_str is None:
            shading_str = False

        with _pygmt.config(FONT_TITLE=titlesize, FONT_LABEL=axes_labelsize,
                           FONT_ANNOT=tick_labelsize):
            _pygmt.makecpt(series=cmap_limits, cmap=cmap, reverse=cmap_reverse,
                           continuous=cmap_continuous)
            figure.shift_origin(xshift=xshift, yshift=yshift)
            figure.grdimage(self.to_xarray(), region=region,
                            projection=proj_str, frame=frame,
                            shading=shading_str)
            if colorbar is not None:
                if shading is not None:
                    figure.colorbar(position=position, frame=cb_str,
                                    shading=shading_amplitude)
                else:
                    figure.colorbar(position=position, frame=cb_str)

        if isinstance(shading, SHGrid):
            f.close()

        return figure


# ---- Complex Driscoll and Healy grid class ----
class DHComplexGrid(SHGrid):
    """
    Class for complex Driscoll and Healy (1994) grids.
    """
    @staticmethod
    def istype(kind):
        return kind == 'complex'

    @staticmethod
    def isgrid(grid):
        return grid == 'DH'

    def __init__(self, array, units=None, copy=True):
        self.nlat, self.nlon = array.shape

        if self.nlat % 2 != 0:
            self.n = self.nlat - 1
            self.extend = True
        else:
            self.n = self.nlat
            self.extend = False

        if self.nlon == 2 * self.nlat - self.extend:
            self.sampling = 2
        elif self.nlon == self.nlat:
            self.sampling = 1
        else:
            raise ValueError('Input array has shape (nlat={:d}, nlon={:d}) '
                             .format(self.nlat, self.nlon) +
                             'but needs nlon=nlat, nlon=2*nlat, or '
                             'nlon=2*nlat-1.'
                             )

        self.lmax = int(self.n / 2 - 1)
        self.grid = 'DH'
        self.kind = 'complex'
        self.units = units

        if copy:
            self.data = _np.copy(array)
        else:
            self.data = array

    def _lats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        if self.extend:
            lats = _np.linspace(90.0, -90.0, num=self.nlat)
        else:
            lats = _np.linspace(90.0, -90.0 + 180.0 / self.nlat, num=self.nlat)
        return lats

    def _lons(self):
        """
        Return a vector containing the longitudes (in degrees) of each row
        of the gridded data.
        """
        if self.extend:
            lons = _np.linspace(0.0, 360.0, num=self.nlon)
        else:
            lons = _np.linspace(0.0, 360.0 - 360.0 / self.nlon, num=self.nlon)
        return lons

    def _histogram(self, **kwargs):
        """Return an area-weighted histogram normalized to unity."""
        raise NotImplementedError('histogram() does not support complex data.')

    def _expand(self, normalization, csphase, lmax_calc, backend, nthreads):
        """Expand the grid into real spherical harmonics."""
        from .shcoeffs import SHCoeffs
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'schmidt':
            norm = 3
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        cilm = backend_module(
            backend=backend, nthreads=nthreads).SHExpandDHC(
                self.data[:self.nlat-self.extend, :self.nlon-self.extend],
                norm=norm, csphase=csphase, sampling=self.sampling,
                lmax_calc=lmax_calc)
        coeffs = SHCoeffs.from_array(cilm, normalization=normalization.lower(),
                                     csphase=csphase, units=self.units,
                                     copy=False)
        return coeffs

    def _plot(self, projection=None, xlabel=None, ylabel=None, colorbar=None,
              cb_triangles=None, cb_label=None, grid=False, ticks=None,
              axes_labelsize=None, tick_labelsize=None, title=None,
              titlesize=None, cmap=None, ax=None, ax2=None,
              tick_interval=None, minor_tick_interval=None, cb_ylabel=None,
              cb_tick_interval=None, cb_minor_tick_interval=None,
              cmap_limits=None, cmap_reverse=None, cmap_limits_complex=None,
              cb_offset=None, cb_width=None):
        """Plot the raw data as a matplotlib simple cylindrical projection,
           or with Cartopy when projection is specified."""
        if ax is None:
            if colorbar is not None:
                if colorbar in set(['top', 'bottom']):
                    scale = 1.5
                else:
                    scale = 1.1
            else:
                scale = 1.2
            figsize = (_mpl.rcParams['figure.figsize'][0],
                       _mpl.rcParams['figure.figsize'][0]*scale)
            fig, axes = _plt.subplots(2, 1, figsize=figsize)
            axreal = axes.flat[0]
            axcomplex = axes.flat[1]
        else:
            axreal = ax
            axcomplex = ax2

        self.to_real().plot(projection=projection, tick_interval=tick_interval,
                            minor_tick_interval=minor_tick_interval,
                            colorbar=colorbar, cb_triangles=cb_triangles,
                            cb_label=cb_label, ticks=ticks,
                            cb_tick_interval=cb_tick_interval,
                            cb_minor_tick_interval=cb_minor_tick_interval,
                            grid=grid, axes_labelsize=axes_labelsize,
                            tick_labelsize=tick_labelsize, cb_offset=cb_offset,
                            title=title[0], titlesize=titlesize,
                            xlabel=xlabel, ylabel=ylabel, cb_ylabel=cb_ylabel,
                            cb_width=cb_width, cmap=cmap,
                            cmap_limits=cmap_limits, cmap_reverse=cmap_reverse,
                            ax=axreal)

        self.to_imag().plot(projection=projection, tick_interval=tick_interval,
                            minor_tick_interval=minor_tick_interval,
                            colorbar=colorbar, cb_triangles=cb_triangles,
                            cb_label=cb_label, ticks=ticks,
                            cb_tick_interval=cb_tick_interval,
                            cb_minor_tick_interval=cb_minor_tick_interval,
                            grid=grid, axes_labelsize=axes_labelsize,
                            tick_labelsize=tick_labelsize, cb_ylabel=cb_ylabel,
                            title=title[1], titlesize=titlesize,
                            cmap=cmap, cmap_limits=cmap_limits_complex,
                            cmap_reverse=cmap_reverse, cb_offset=cb_offset,
                            cb_width=cb_width, xlabel=xlabel, ylabel=ylabel,
                            ax=axcomplex)

        if ax is None:
            return fig, axes

    def _plot_pygmt(self, fig=None, projection=None, region=None, width=None,
                    unit=None, central_latitude=None, central_longitude=None,
                    grid=None, tick_interval=None, minor_tick_interval=None,
                    ticks=None, title=None, cmap=None, cmap_limits=None,
                    cmap_limits_complex=None, cmap_reverse=None,
                    cmap_continuous=None, colorbar=None, cb_triangles=None,
                    cb_label=None, cb_ylabel=None, cb_tick_interval=None,
                    cb_minor_tick_interval=None, titlesize=None,
                    axes_labelsize=None, tick_labelsize=None, horizon=None,
                    offset=None, cb_offset=None):
        """
        Plot projected data using pygmt.
        """
        if fig is None:
            figure = _pygmt.Figure()
        else:
            figure = fig

        self.to_imag().plotgmt(fig=figure, projection=projection,
                               region=region, width=width, unit=unit,
                               central_latitude=central_latitude,
                               central_longitude=central_longitude,
                               grid=grid, tick_interval=tick_interval,
                               minor_tick_interval=minor_tick_interval,
                               ticks=ticks, title=title[1], cmap=cmap,
                               cmap_limits=cmap_limits_complex,
                               cmap_reverse=cmap_reverse, cb_offset=cb_offset,
                               cmap_continuous=cmap_continuous,
                               colorbar=colorbar, cb_triangles=cb_triangles,
                               cb_label=cb_label, cb_ylabel=cb_ylabel,
                               cb_tick_interval=cb_tick_interval,
                               cb_minor_tick_interval=cb_minor_tick_interval,
                               titlesize=titlesize,
                               axes_labelsize=axes_labelsize,
                               tick_labelsize=tick_labelsize, horizon=horizon,
                               offset=offset)

        offset_real = _np.copy(offset)
        if offset_real[1] is None:
            offset_real[1] = width / 2. + 50. / 72.
        else:
            offset_real[1] += width / 2. + 50. / 72.

        self.to_real().plotgmt(fig=figure, projection=projection,
                               region=region, width=width, unit=unit,
                               central_latitude=central_latitude,
                               central_longitude=central_longitude,
                               grid=grid, tick_interval=tick_interval,
                               minor_tick_interval=minor_tick_interval,
                               ticks=ticks, title=title[0], cmap=cmap,
                               cmap_limits=cmap_limits, cb_offset=cb_offset,
                               cmap_reverse=cmap_reverse,
                               cmap_continuous=cmap_continuous,
                               colorbar=colorbar, cb_triangles=cb_triangles,
                               cb_label=cb_label, cb_ylabel=cb_ylabel,
                               cb_tick_interval=cb_tick_interval,
                               cb_minor_tick_interval=cb_minor_tick_interval,
                               titlesize=titlesize,
                               axes_labelsize=axes_labelsize,
                               tick_labelsize=tick_labelsize,
                               horizon=horizon, offset=offset_real)

        return figure


# ---- Real Gauss-Legendre Quadrature grid class ----

class GLQRealGrid(SHGrid):
    """
    Class for real Gauss-Legendre Quadrature grids.
    """
    @staticmethod
    def istype(kind):
        return kind == 'real'

    @staticmethod
    def isgrid(grid):
        return grid == 'GLQ'

    def __init__(self, array, zeros=None, weights=None, units=None, copy=True):
        self.nlat, self.nlon = array.shape
        self.lmax = self.nlat - 1

        if self.nlon == 2 * self.lmax + 1:
            self.extend = False
        elif self.nlon == 2 * self.lmax + 2:
            self.extend = True
        else:
            raise ValueError('Input array has shape (nlat={:d}, nlon={:d}) '
                             .format(self.nlat, self.nlon) +
                             'but needs ({:d}, {:d}) or ({:d}, {:d}).'
                             .format(self.lmax+1, 2*self.lmax+1, self.lmax+1,
                                     2*self.lmax+2)
                             )

        if zeros is None or weights is None:
            self.zeros, self.weights = _shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights = weights

        self.grid = 'GLQ'
        self.kind = 'real'
        self.units = units

        if copy:
            self.data = _np.copy(array)
        else:
            self.data = array

    def _lats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = 90. - _np.arccos(self.zeros) * 180. / _np.pi
        return lats

    def _lons(self):
        """
        Return a vector containing the longitudes (in degrees) of each column
        of the gridded data.
        """
        if self.extend:
            lons = _np.linspace(0.0, 360.0, num=self.nlon)
        else:
            lons = _np.linspace(0.0, 360.0 - 360.0 / self.nlon, num=self.nlon)
        return lons

    def _histogram(self, bins=None, range=None):
        """Return an area-weighted histogram normalized to unity."""
        delta_phi = self.lons()[1] - self.lons()[0]
        delta_phi *= _np.pi / 180.

        da = _np.zeros_like(self.data)

        # First latitudinal band includes 90 N
        theta2 = 90.0 - (self.lats()[0] + self.lats()[1]) / 2.
        da[0, :] = 1. - _np.cos(theta2 * _np.pi / 180.)
        for i in _np.arange(1, self.nlat-1):
            theta1 = 90.0 - (self.lats()[i-1] + self.lats()[i]) / 2.
            theta2 = 90.0 - (self.lats()[i] + self.lats()[i+1]) / 2.
            da[i, :] = _np.cos(theta1 * _np.pi / 180.) - \
                _np.cos(theta2 * _np.pi / 180.)
        # Last latitudinal band includes 90 S
        i = self.nlat - 1
        theta1 = 90.0 - (self.lats()[i] + self.lats()[i-1]) / 2.
        da[i, :] = _np.cos(theta1 * _np.pi / 180.) + 1

        return _np.histogram(self.data[:, :self.nlon-self.extend],
                             bins=bins,
                             weights=da[:, :self.nlon-self.extend],
                             density=True, range=range)

    def _expand(self, normalization, csphase, lmax_calc, backend, nthreads):
        """Expand the grid into real spherical harmonics."""
        from .shcoeffs import SHCoeffs
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'unnorm':
            norm = 3
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt' " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        cilm = backend_module(
            backend=backend, nthreads=nthreads).SHExpandGLQ(
                self.data[:, :self.nlon-self.extend], self.weights, self.zeros,
                norm=norm, csphase=csphase, lmax_calc=lmax_calc)
        coeffs = SHCoeffs.from_array(cilm, normalization=normalization.lower(),
                                     csphase=csphase, units=self.units,
                                     copy=False)
        return coeffs

    def _plot(self, projection=None, xlabel=None, ylabel=None, ax=None,
              ax2=None, colorbar=None, cb_triangles=None, cb_label=None,
              grid=False, axes_labelsize=None, tick_labelsize=None,
              title=None, titlesize=None, cmap=None, tick_interval=None,
              minor_tick_interval=None, cb_tick_interval=None, ticks=None,
              cb_minor_tick_interval=None, cmap_limits=None, cmap_reverse=None,
              cmap_limits_complex=None, cb_ylabel=None, cb_offset=None,
              cb_width=None):
        """Plot the data using a matplotlib cylindrical projection."""
        if ax is None:
            if colorbar is not None:
                if colorbar in set(['top', 'bottom']):
                    scale = 0.67
                else:
                    scale = 0.5
            else:
                scale = 0.55
            figsize = (_mpl.rcParams['figure.figsize'][0],
                       _mpl.rcParams['figure.figsize'][0] * scale)
            fig, axes = _plt.subplots(1, 1, figsize=figsize)
        else:
            axes = ax

        if tick_interval[0] is None:
            xticks = []
        else:
            xticks = _np.arange(0, self.nlon, tick_interval[0])
        if tick_interval[1] is None:
            yticks = []
        else:
            yticks = _np.arange(0, self.nlat, tick_interval[1])
        if minor_tick_interval[0] is None:
            minor_xticks = []
        else:
            minor_xticks = _np.arange(0, self.nlon, minor_tick_interval[0])
        if minor_tick_interval[1] is None:
            minor_yticks = []
        else:
            minor_yticks = _np.arange(0, self.nlat, minor_tick_interval[1])

        # make colormap
        if cmap_limits is None:
            cmap_limits = [self.min(), self.max()]
        if len(cmap_limits) == 3:
            num = int((cmap_limits[1] - cmap_limits[0]) / cmap_limits[2])
            if isinstance(cmap, _mpl.colors.Colormap):
                cmap_scaled = cmap._resample(num)
            else:
                cmap_scaled = _mpl.cm.get_cmap(cmap, num)
        else:
            cmap_scaled = _mpl.cm.get_cmap(cmap)
        if cmap_reverse:
            cmap_scaled = cmap_scaled.reversed()

        # compute colorbar ticks
        cb_ticks = None
        cb_minor_ticks = None
        vmin = cmap_limits[0]
        vmax = cmap_limits[1]
        if cb_tick_interval is not None:
            if _np.sign(vmin) == -1.:
                start = (abs(vmin) // cb_tick_interval) \
                    * cb_tick_interval * _np.sign(vmin)
            else:
                start = (vmin // cb_tick_interval + 1) \
                    * cb_tick_interval
            if _np.sign(vmax) == -1.:
                stop = (abs(vmax) // cb_tick_interval + 1) \
                    * cb_tick_interval * _np.sign(vmax)
            else:
                stop = (vmax // cb_tick_interval) * cb_tick_interval
            cb_ticks = _np.linspace(start, stop,
                                    num=int((stop-start)//cb_tick_interval+1),
                                    endpoint=True)
        if cb_minor_tick_interval is not None:
            if _np.sign(vmin) == -1.:
                start = (abs(vmin) // cb_minor_tick_interval) \
                    * cb_minor_tick_interval * _np.sign(vmin)
            else:
                start = (vmin // cb_minor_tick_interval + 1) \
                    * cb_minor_tick_interval
            if _np.sign(vmax) == -1.:
                stop = (abs(vmax) // cb_minor_tick_interval + 1) \
                    * cb_minor_tick_interval * _np.sign(vmax)
            else:
                stop = (vmax // cb_minor_tick_interval) * \
                    cb_minor_tick_interval
            cb_minor_ticks = _np.linspace(
                start, stop, num=int((stop-start)//cb_minor_tick_interval+1),
                endpoint=True)

        # determine which ticks to plot
        if 'W' in ticks:
            left, labelleft = True, True
        elif 'w' in ticks:
            left, labelleft = True, False
        else:
            left, labelleft = False, False
        if 'S' in ticks:
            bottom, labelbottom = True, True
        elif 's' in ticks:
            bottom, labelbottom = True, False
        else:
            bottom, labelbottom = False, False
        if 'E' in ticks:
            right, labelright = True, True
        elif 'e' in ticks:
            right, labelright = True, False
        else:
            right, labelright = False, False
        if 'N' in ticks:
            top, labeltop = True, True
        elif 'n' in ticks:
            top, labeltop = True, False
        else:
            top, labeltop = False, False

        # plot image, ticks, and annotations
        extent = (-0.5, self.nlon-0.5, -0.5, self.nlat-0.5)
        cim = axes.imshow(self.data, extent=extent, origin='upper',
                          cmap=cmap_scaled, vmin=cmap_limits[0],
                          vmax=cmap_limits[1])
        axes.set(xticks=xticks, yticks=yticks)
        axes.set_xlabel(xlabel, fontsize=axes_labelsize)
        axes.set_ylabel(ylabel, fontsize=axes_labelsize)
        axes.set_xticklabels(xticks, fontsize=tick_labelsize)
        axes.set_yticklabels(yticks, fontsize=tick_labelsize)
        axes.set_xticks(minor_xticks, minor=True)
        axes.set_yticks(minor_yticks, minor=True)
        axes.tick_params(bottom=bottom, top=top, right=right, left=left,
                         labelbottom=labelbottom, labeltop=labeltop,
                         labelleft=labelleft, labelright=labelright,
                         which='both')
        axes.grid(grid, which='major')
        if title is not None:
            axes.set_title(title, fontsize=titlesize)

        # plot colorbar
        if colorbar is not None:
            if cb_offset is None:
                if colorbar in set(['left', 'right']):
                    offset = 0.15
                    if (colorbar == 'left' and 'W' in ticks) or \
                            (colorbar == 'right' and 'E' in ticks):
                        offset += 2 * tick_labelsize / 72.
                    # add space for ylabel on left of plot only
                    if ylabel != '' and ylabel is not None and \
                            colorbar == 'left':
                        offset += 1.4 * axes_labelsize / 72.
                else:
                    offset = 0.
                    # add space for ticks
                    if (colorbar == 'bottom' and bottom) or \
                            (colorbar == 'top' and top):
                        offset += _mpl.rcParams['xtick.major.size']
                    # add space for labels
                    if (colorbar == 'bottom' and labelbottom) or \
                            (colorbar == 'top' and labeltop):
                        offset += _mpl.rcParams['xtick.major.pad']
                        offset += tick_labelsize
                    # add space for xlabel on bottom of plot only
                    if xlabel != '' and xlabel is not None and \
                            colorbar == 'bottom':
                        offset += axes_labelsize
                    offset += 1.3 * _mpl.rcParams['font.size']  # add extra
                    offset /= 72.  # convert to inches
            else:
                offset = cb_offset / 72.  # convert to inches

            divider = _make_axes_locatable(axes)
            if colorbar in set(['left', 'right']):
                orientation = 'vertical'
                if cb_width is None:
                    size = '2.5%'
                else:
                    size = '{:f}%'.format(cb_width)
            else:
                orientation = 'horizontal'
                if cb_width is None:
                    size = '5%'
                else:
                    size = '{:f}%'.format(cb_width)
            cax = divider.append_axes(colorbar, size=size, pad=offset,
                                      axes_class=_plt.Axes)
            cbar = _plt.colorbar(cim, cax=cax, orientation=orientation,
                                 extend=cb_triangles)
            if colorbar == 'left':
                cbar.ax.yaxis.set_ticks_position('left')
                cbar.ax.yaxis.set_label_position('left')
            if colorbar == 'top':
                cbar.ax.xaxis.set_ticks_position('top')
                cbar.ax.xaxis.set_label_position('top')
            if cb_label is not None:
                cbar.set_label(cb_label, fontsize=axes_labelsize)
            if cb_ylabel is not None:
                if colorbar in set(['left', 'right']):
                    cbar.ax.xaxis.set_label_position('top')
                    cbar.ax.set_xlabel(cb_ylabel, fontsize=tick_labelsize)
                else:
                    cbar.ax.yaxis.set_label_position('right')
                    cbar.ax.set_ylabel(cb_ylabel, fontsize=tick_labelsize,
                                       rotation=0., labelpad=axes_labelsize/2.,
                                       va='center', ha='left')
            if cb_ticks is not None:
                cbar.set_ticks(cb_ticks)
            if cb_minor_ticks is not None:
                if colorbar in set(['top', 'bottom']):
                    cbar.ax.xaxis.set_ticks(cb_minor_ticks, minor=True)
                else:
                    cbar.ax.yaxis.set_ticks(cb_minor_ticks, minor=True)
            cbar.ax.tick_params(labelsize=tick_labelsize)

        if ax is None:
            return fig, axes

    def _plot_pygmt(self, **kwargs):
        """
        Plot the projected data using pygmt.
        """
        raise NotImplementedError('plotgmt() does not support the plotting '
                                  'of GLQ gridded data.')


# ---- Complex Gauss-Legendre Quadrature grid class ----

class GLQComplexGrid(SHGrid):
    """
    Class for complex Gauss-Legendre Quadrature grids.
    """
    @staticmethod
    def istype(kind):
        return kind == 'complex'

    @staticmethod
    def isgrid(grid):
        return grid == 'GLQ'

    def __init__(self, array, zeros=None, weights=None, units=None, copy=True):
        self.nlat, self.nlon = array.shape
        self.lmax = self.nlat - 1

        if self.nlon == 2 * self.lmax + 1:
            self.extend = False
        elif self.nlon == 2 * self.lmax + 2:
            self.extend = True
        else:
            raise ValueError('Input array has shape (nlat={:d}, nlon={:d}) '
                             .format(self.nlat, self.nlon) +
                             'but needs ({:d}, {:d}) or ({:d}, {:d}).'
                             .format(self.lmax+1, 2*self.lmax+1, self.lmax+1,
                                     2*self.lmax+2)
                             )

        if zeros is None or weights is None:
            self.zeros, self.weights = _shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights = weights

        self.grid = 'GLQ'
        self.kind = 'complex'
        self.units = units

        if copy:
            self.data = _np.copy(array)
        else:
            self.data = array

    def _lats(self):
        """Return the latitudes (in degrees) of the gridded data rows."""
        lats = 90. - _np.arccos(self.zeros) * 180. / _np.pi
        return lats

    def _lons(self):
        """Return the longitudes (in degrees) of the gridded data columns."""
        if self.extend:
            lons = _np.linspace(0.0, 360.0, num=self.nlon)
        else:
            lons = _np.linspace(0.0, 360.0 - 360.0 / self.nlon, num=self.nlon)
        return lons

    def _histogram(self, **kwargs):
        """Return an area-weighted histogram normalized to unity."""
        raise NotImplementedError('histogram() does not support complex data.')

    def _expand(self, normalization, csphase, lmax_calc, backend, nthreads):
        """Expand the grid into real spherical harmonics."""
        from .shcoeffs import SHCoeffs
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'unnorm':
            norm = 3
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt' " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        cilm = backend_module(
            backend=backend, nthreads=nthreads).SHExpandGLQC(
                self.data[:, :self.nlon-self.extend], self.weights, self.zeros,
                norm=norm, csphase=csphase, lmax_calc=lmax_calc)
        coeffs = SHCoeffs.from_array(cilm, normalization=normalization.lower(),
                                     csphase=csphase, units=self.units,
                                     copy=False)
        return coeffs

    def _plot(self, projection=None, minor_xticks=[], minor_yticks=[],
              xlabel=None, ylabel=None, ax=None, ax2=None, colorbar=None,
              cb_triangles=None, cb_label=None, grid=False, ticks=None,
              axes_labelsize=None, tick_labelsize=None, title=None,
              titlesize=None, cmap=None, tick_interval=None, cb_ylabel=None,
              minor_tick_interval=None, cb_tick_interval=None,
              cb_minor_tick_interval=None, cmap_limits=None, cmap_reverse=None,
              cmap_limits_complex=None, cb_offset=None, cb_width=None):
        """Plot the raw data using a simply cylindrical projection."""
        if ax is None:
            if colorbar is not None:
                if colorbar in set(['top', 'bottom']):
                    scale = 1.5
                else:
                    scale = 1.1
            else:
                scale = 1.2
            figsize = (_mpl.rcParams['figure.figsize'][0],
                       _mpl.rcParams['figure.figsize'][0]*scale)
            fig, axes = _plt.subplots(2, 1, figsize=figsize)
            axreal = axes.flat[0]
            axcomplex = axes.flat[1]
        else:
            axreal = ax
            axcomplex = ax2

        self.to_real().plot(projection=projection, tick_interval=tick_interval,
                            minor_tick_interval=minor_tick_interval,
                            colorbar=colorbar, cb_triangles=cb_triangles,
                            cb_label=cb_label, ticks=ticks,
                            cb_tick_interval=cb_tick_interval,
                            cb_minor_tick_interval=cb_minor_tick_interval,
                            grid=grid, axes_labelsize=axes_labelsize,
                            tick_labelsize=tick_labelsize, cb_offset=cb_offset,
                            title=title[0], titlesize=titlesize,
                            xlabel=xlabel, ylabel=ylabel, cb_ylabel=cb_ylabel,
                            cb_width=cb_width, cmap=cmap,
                            cmap_limits=cmap_limits, cmap_reverse=cmap_reverse,
                            ax=axreal)

        self.to_imag().plot(projection=projection, tick_interval=tick_interval,
                            minor_tick_interval=minor_tick_interval,
                            colorbar=colorbar, cb_triangles=cb_triangles,
                            cb_label=cb_label, ticks=ticks,
                            cb_tick_interval=cb_tick_interval,
                            cb_minor_tick_interval=cb_minor_tick_interval,
                            grid=grid, axes_labelsize=axes_labelsize,
                            tick_labelsize=tick_labelsize, cb_offset=cb_offset,
                            title=title[1], titlesize=titlesize,
                            cmap=cmap, cmap_limits=cmap_limits_complex,
                            cmap_reverse=cmap_reverse, cb_ylabel=cb_ylabel,
                            cb_width=cb_width, xlabel=xlabel, ylabel=ylabel,
                            ax=axcomplex)

        if ax is None:
            return fig, axes

    def _plot_pygmt(self, **kwargs):
        """
        Plot the projected data using pygmt.
        """
        raise NotImplementedError('plotgmt() does not support the plotting '
                                  'of complex GLQ gridded data.')
