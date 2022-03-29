"""
    Class for the height of the geoid.
"""
import numpy as _np
import copy as _copy
import xarray as _xr

from .shgrid import SHGrid as _SHGrid
from .shgrid import _pygmt_module


class SHGeoid(object):
    """
    Class for the height of the geoid. The class is initialized from a class
    instance of SHGravCoeffs using the method geoid(). Geoid heights are
    referenced to a flattened ellipsoid of semimajor axis a and flattening f.

    Attributes:

    geoid          : SHGrid class instance of the geoid.
    gm             : Gravitational constant time the mass of the body.
    potref         : Potential of the chosen geoid, in m2/s2.
    a              : Semimajor axis of the reference ellipsoid.
    f              : Flattening of the reference ellipsoid, f=(a-b)/a.
    omega          : Angular rotation rate of the body.
    r              : Reference radius of the Taylor expansion.
    order          : Order of the Taylor expansion.
    lmax           : The maximum spherical harmonic degree resolvable by the
                     geoid grid.
    lmax_calc      : The maximum spherical harmonic degree of the gravitational
                     potential used in creating the geoid.
    nlat, nlon     : The number of latitude and longitude bands in the geoid
                     grid.
    n              : The number of samples in latitude.
    sampling       : The longitudinal sampling for Driscoll and Healy grids.
                     Either 1 for equally sampled grids (nlat=nlon) or 2 for
                     equally spaced grids in degrees.
    units          : The units of the gridded data.
    extend         : True if the grid contains the redundant column for 360 E
                     and the unnecessary row for 90 S.
    epoch          : The epoch time of the gravity model.

    Methods:

    plot()         : Plot the geoid.
    to_xarray()    : Return the gridded data as an xarray DataArray.
    to_netcdf()    : Return the gridded data as a netcdf formatted file or
                     object.
    copy()         : Return a copy of the class instance.
    info()         : Print a summary of the data stored in the SHGrid instance.
    """
    def __init__(self, geoid, gm, potref, a, f, omega, r, order, lmax,
                 lmax_calc, units=None, epoch=None):
        """
        Initialize the SHGeoid class.
        """
        self.geoid = _SHGrid.from_array(geoid, grid='DH', units=units)
        self.grid = self.geoid.grid
        self.sampling = self.geoid.sampling
        self.nlat = self.geoid.nlat
        self.nlon = self.geoid.nlon
        self.n = self.geoid.n
        self.extend = self.geoid.extend
        self.potref = potref
        self.gm = gm
        self.a = a
        self.f = f
        if omega is None:
            self.omega = 0.0
        else:
            self.omega = omega
        self.order = order
        self.r = r
        self.lmax = lmax
        self.lmax_calc = lmax_calc
        self.units = units
        self.epoch = epoch

    def copy(self):
        """
        Return a deep copy of the class instance.

        Usage
        -----
        copy = x.copy()
        """
        return _copy.deepcopy(self)

    def info(self):
        """
        Print a summary of the data stored in the SHGeoid class instance.

        Usage
        -----
        x.info()
        """
        print(repr(self))

    def __repr__(self):
        str = ('grid = {:s}\n'
               'sampling = {:d}\n'
               'nlat = {:d}\n'
               'nlon = {:d}\n'
               'n = {:d}\n'
               'sampling = {:d}\n'
               'extend = {}\n'
               'lmax = {:d}\n'
               'lmax_calc = {:d}\n'
               'gm (m3 / s2) = {:e}\n'
               'reference potential (m2 /s2) = {:e}\n'
               'a (m)= {:e}\n'
               'f = {:e}\n'
               'omega (rad / s) = {:s}\n'
               'radius of Taylor expansion (m) = {:e}\n'
               'order of expansion = {:d}\n'
               'units = {:s}\n'
               'epoch = {:s}'
               .format(repr(self.grid), self.sampling, self.nlat, self.nlon,
                       self.n, self.sampling, self.extend, self.lmax,
                       self.lmax_calc, self.gm, self.potref, self.a, self.f,
                       repr(self.omega), self.r, self.order, repr(self.units),
                       repr(self.epoch)))
        return str

    def plot(self, projection=None, tick_interval=[30, 30],
             minor_tick_interval=[None, None], xlabel=None, ylabel=None,
             title=None, titlesize=None, colorbar='right', cmap='viridis',
             cmap_limits=None, cmap_reverse=False, cb_triangles='neither',
             cb_label='geoid, m', cb_tick_interval=None, grid=False,
             axes_labelsize=None, tick_labelsize=None, show=True, ax=None,
             cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
             fname=None, cb_offset=None, cb_width=None):
        """
        Plot the geoid.

        Usage
        -----
        x.plot([projection, tick_interval, minor_tick_interval, ticks,
                xlabel, ylabel, title, colorbar, cmap, cmap_limits,
                cmap_reverse, cb_triangles, cb_label, cb_ylabel,
                cb_tick_interval, cb_minor_tick_interval, cb_offset, cb_width,
                grid, titlesize, axes_labelsize, tick_labelsize, ax, show,
                fname])

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
            and north boundaries of the plot.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        title : str or list, optional, default = None
            The title of the plot.
        colorbar : str, optional, default = 'right'
            Plot a colorbar along the 'top', 'right', 'bottom', or 'left' axis.
        cmap : str, optional, default = 'viridis'
            The color map to use when plotting the data and colorbar.
        cmap_limits : list, optional, default = [self.min(), self.max()]
            Set the lower and upper limits of the data used by the colormap,
            and optionally an interval for each color band. If the interval is
            specified, the number of discrete colors will be
            (cmap_limits[1]-cmap_limits[0])/cmap_limits[2].
        cmap_reverse : bool, optional, default = False
            Set to True to reverse the sense of the color progression in the
            color table.
        cb_triangles : str, optional, default = 'neither'
            Add triangles to the edges of the colorbar for minimum and maximum
            values. Can be 'neither', 'both', 'min', or 'max'.
        cb_label : str, optional, default = 'geoid, m'
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
        titlesize : int, optional, default = None
            The font size of the title.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
        """
        return self.geoid.plot(projection=projection,
                               tick_interval=tick_interval,
                               minor_tick_interval=minor_tick_interval,
                               xlabel=xlabel, ylabel=ylabel, title=title,
                               titlesize=titlesize, colorbar=colorbar,
                               cmap=cmap, cmap_limits=cmap_limits,
                               cmap_reverse=cmap_reverse, cb_offset=cb_offset,
                               cb_triangles=cb_triangles, cb_label=cb_label,
                               cb_tick_interval=cb_tick_interval, grid=grid,
                               cb_ylabel=cb_ylabel, ticks=ticks,
                               cb_minor_tick_interval=cb_minor_tick_interval,
                               axes_labelsize=axes_labelsize,
                               cb_width=cb_width,
                               tick_labelsize=tick_labelsize, ax=ax,
                               show=show, fname=fname)

    def to_netcdf(self, filename=None, title='', description='',
                  comment='pyshtools grid', name='geoid',
                  dtype='d'):
        """
        Return the gridded data as a netcdf formatted file or object.

        Usage
        -----
        x.to_netcdf([filename, title, description, comment, name, dtype])

        Parameters
        ----------
        filename : str, optional, default = None
            Name of output file.
        title : str, optional, default = ''
            Title of the dataset.
        description : str, optional, default = ''
            Description of the dataset ('Remark' in gmt grd files).
        comment : str, optional, default = 'pyshtools grid'
            Additional information about how the data were generated.
        name : str, optional, default = 'data'
            Name of the data array.
        dtype : str, optional, default = 'd'
            Data type of the output array. Either 'f' or 'd' for single or
            double precision floating point, respectively.
        """
        if dtype == 'f':
            _nparray = self.geoid.to_array().astype(_np.float32)
        elif dtype == 'd':
            _nparray = self.geoid.to_array()
        else:
            raise ValueError("dtype must be either 'f' or 'd' for single or "
                             "double precision floating point.")

        attrs = {'actual_range': [self.geoid.min(), self.geoid.max()],
                 'title': title,
                 'comment': comment,
                 'long_name': 'Geoid',
                 'units': self.units,
                 'nlat': self.geoid.nlat,
                 'nlon': self.geoid.nlon,
                 'lmax': self.lmax,
                 'kind': self.geoid.kind,
                 'grid': self.geoid.grid,
                 'gm': self.gm,
                 'potref': self.potref,
                 'a': self.a,
                 'f': self.f,
                 'omega': self.omega,
                 'r': self.r,
                 'order': self.order,
                 'lmax_calc': self.lmax_calc,
                 'sampling': self.geoid.sampling,
                 'n': self.n,
                 'extend': repr(self.extend)
                 }
        if self.epoch is not None:
            attrs['epoch'] = self.epoch

        _data = _xr.DataArray(_nparray, dims=('latitude', 'longitude'),
                              coords=[('latitude', self.geoid.lats(),
                                       {'units': 'degrees_north'}),
                                      ('longitude', self.geoid.lons(),
                                       {'units': 'degrees_east'})],
                              attrs=attrs)
        _dataset = _xr.Dataset({name: _data},
                               attrs={'title': title,
                                      'comment': comment,
                                      'description': description})
        if filename is None:
            return _dataset.to_netcdf()
        else:
            _dataset.to_netcdf(filename)

    def to_xarray(self, title='', comment='pyshtools grid'):
        """
        Return the gridded data as an xarray DataArray.

        Usage
        -----
        x.to_xarray([title, comment])

        Parameters
        ----------
        title : str, optional, default = ''
            Title of the dataset.
        comment : str, optional, default = 'pyshtools grid'
            Additional information about how the data were generated.
        """
        attrs = {'actual_range': [self.geoid.min(), self.geoid.max()],
                 'title': title,
                 'comment': comment,
                 'long_name': 'Geoid',
                 'units': self.units,
                 'nlat': self.geoid.nlat,
                 'nlon': self.geoid.nlon,
                 'lmax': self.lmax,
                 'kind': self.geoid.kind,
                 'grid': self.geoid.grid,
                 'gm': self.gm,
                 'potref': self.potref,
                 'a': self.a,
                 'f': self.f,
                 'omega': self.omega,
                 'r': self.r,
                 'order': self.order,
                 'lmax_calc': self.lmax_calc,
                 'sampling': self.geoid.sampling,
                 'n': self.n,
                 'extend': repr(self.extend)
                 }
        if self.epoch is not None:
            attrs['epoch'] = self.epoch

        da = _xr.DataArray(self.geoid.to_array(),
                           dims=('latitude', 'longitude'),
                           coords=[('latitude', self.geoid.lats(),
                                    {'units': 'degrees_north'}),
                                   ('longitude', self.geoid.lons(),
                                    {'units': 'degrees_east'})],
                           attrs=attrs)
        if _pygmt_module:
            da.gmt.registration = 0
            da.gmt.gtype = 1
        return da
