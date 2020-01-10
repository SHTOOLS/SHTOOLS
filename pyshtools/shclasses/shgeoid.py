"""
    Class for the height of the geoid.
"""
import numpy as _np
import copy as _copy
import xarray as _xr

from .shcoeffsgrid import SHGrid as _SHGrid


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
    extend         : True if the grid contains the redundant column for 360 E
                     and the unnecessary row for 90 S.

    Methods:

    plot()         : Plot the geoid.
    to_xarray()    : Return the gridded data as an xarray DataArray.
    to_netcdf()    : Return the gridded data as a netcdf formatted file or
                     object.
    copy()         : Return a copy of the class instance.
    info()         : Print a summary of the data stored in the SHGrid instance.
    """
    def __init__(self, geoid, gm, potref, a, f, omega, r, order, lmax,
                 lmax_calc):
        """
        Initialize the SHGeoid class.
        """
        self.geoid = _SHGrid.from_array(geoid, grid='DH')
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
        str = ('grid = {:s}\n'.format(repr(self.grid)))
        if self.grid == 'DH':
            str += 'sampling = {:d}\n'.format(self.sampling)
        str += ('nlat = {:d}\n'
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
                'order of expansion = {:d}'
                .format(self.nlat, self.nlon, self.n, self.sampling,
                        self.extend, self.lmax, self.lmax_calc, self.gm,
                        self.potref, self.a, self.f, repr(self.omega), self.r,
                        self.order))
        return str

    def plot(self, tick_interval=[30, 30], minor_tick_interval=[None, None],
             colorbar=True, cb_orientation='vertical', cb_label='geoid, m',
             grid=False, axes_labelsize=None, tick_labelsize=None, show=True,
             **kwargs):
        """
        Plot the geoid.

        Usage
        -----
        x.plot([tick_interval, minor_tick_interval, xlabel, ylabel, colorbar,
                cb_orientation, cb_label, grid, axes_labelsize, tick_labelsize,
                ax, show, fname, **kwargs])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = 'geoid, m'
            Text label for the colorbar.
        grid : bool, optional, default = False
            If True, plot major grid lines.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        kwargs : optional
            Keyword arguements that will be sent to plt.imshow(), such as cmap,
            vmin and vmax.
        """
        return self.geoid.plot(tick_interval=tick_interval,
                               minor_tick_interval=minor_tick_interval,
                               colorbar=colorbar,
                               cb_orientation=cb_orientation,
                               cb_label=cb_label,
                               grid=grid, axes_labelsize=axes_labelsize,
                               tick_labelsize=tick_labelsize,
                               show=show, **kwargs)

    def to_netcdf(self, filename=None, title='', description='',
                  comment='Grid generated by pyshtools', name='geoid',
                  dtype='f'):
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
        comment : str, optional, default = 'Grid generated by pyshtools'
            Additional information about how the data were generated.
        name : str, optional, default = 'data'
            Name of the data array.
        dtype : str, optional, default = 'f'
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
                 'units': 'meters',
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
                 'extend': self.extend
                 }

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

    def to_xarray(self, title='', comment='Grid generated by pyshtools'):
        """
        Return the gridded data as an xarray DataArray.

        Usage
        -----
        x.to_xarray([title, comment])

        Parameters
        ----------
        title : str, optional, default = ''
            Title of the dataset.
        comment : str, optional, default = 'Grid generated by pyshtools'
            Additional information about how the data were generated.
        """
        attrs = {'actual_range': [self.geoid.min(), self.geoid.max()],
                 'title': title,
                 'comment': comment,
                 'long_name': 'Geoid',
                 'units': 'meters',
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
                 'extend': self.extend
                 }

        return _xr.DataArray(self.geoid.to_array(),
                             dims=('latitude', 'longitude'),
                             coords=[('latitude', self.geoid.lats(),
                                      {'units': 'degrees_north'}),
                                     ('longitude', self.geoid.lons(),
                                      {'units': 'degrees_east'})],
                             attrs=attrs)
