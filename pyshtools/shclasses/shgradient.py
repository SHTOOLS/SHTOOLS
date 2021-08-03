"""
    Class for grids of the two components of the horizontal gradient.
"""
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy
import xarray as _xr

from .shgrid import SHGrid as _SHGrid


class SHGradient(object):
    """
    Class for grids of the two components of the horizontal gradient of a
    scalar function. The class is initialized from a class instance of
    SHCoeffs using the method gradient().

    Attributes:

    theta          : SHGrid class instance of the theta component of the
                     horizontal gradient.
    phi            : SHGrid class instance of the phi component of the
                     horizontal gradient.
    units          : The units of the gridded data.
    lmax           : The maximum spherical harmonic degree resolvable by the
                     grids.
    lmax_calc      : The maximum spherical harmonic degree of the function
                     used in creating the grids.
    nlat, nlon     : The number of latitude and longitude bands in the grids.
    n              : The number of samples in latitude.
    sampling       : The longitudinal sampling for Driscoll and Healy grids.
                     Either 1 for equally sampled grids (nlat=nlon) or 2 for
                     equally spaced grids in degrees.
    extend         : True if the grid contains the redundant column for 360 E
                     and the unnecessary row for 90 S.

    Methods:

    plot()        : Plot the two components of the horizontal gradient.
    plot_theta()  : Plot the theta component of the horizontal gradient.
    plot_phi()    : Plot the phi component of the horizontal gradient.
    to_xarray()   : Return an xarray DataSet of all gridded data.
    copy()        : Return a copy of the class instance.
    info()        : Print a summary of the data stored in the SHGravGrid
                    instance.
    """

    def __init__(self, theta, phi, lmax, lmax_calc, units=None,
                 pot_units=None, epoch=None):
        """
        Initialize the SHGradient class.
        """
        self.theta = _SHGrid.from_array(theta, grid='DH', units=units)
        self.phi = _SHGrid.from_array(phi, grid='DH', units=units)
        self.grid = self.theta.grid
        self.sampling = self.theta.sampling
        self.nlat = self.theta.nlat
        self.nlon = self.theta.nlon
        self.n = self.theta.n
        self.extend = self.theta.extend
        self.lmax = lmax
        self.lmax_calc = lmax_calc
        self.units = units

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
        Print a summary of the data stored in the SHGradient class instance.

        Usage
        -----
        x.info()
        """
        print(repr(self))

    def __repr__(self):
        str = ('grid = {:s}\n'
               'nlat = {:d}\n'
               'nlon = {:d}\n'
               'n = {:d}\n'
               'sampling = {:d}\n'
               'extend = {}\n'
               'lmax = {:d}\n'
               'lmax_calc = {:d}\n'
               'units = {:s}'
               .format(self.grid, self.nlat, self.nlon, self.n, self.sampling,
                       self.extend, self.lmax, self.lmax_calc,
                       repr(self.units)))
        return str

    def plot_theta(self, projection=None, tick_interval=[30, 30],
                   minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                   title=None, titlesize=None, colorbar='right',
                   cmap='viridis', cmap_limits=None, cmap_reverse=False,
                   cb_triangles='neither', cb_label='$\\theta$ component',
                   cb_tick_interval=None, grid=False, axes_labelsize=None,
                   tick_labelsize=None, show=True, ax=None, cb_offset=None,
                   cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                   fname=None, cb_width=None):
        """
        Plot the theta component of the horizontal gradient.

        Usage
        -----
        x.plot_theta([projection, tick_interval, minor_tick_interval, ticks,
                      xlabel, ylabel, title, colorbar, cmap, cmap_limits,
                      cmap_reverse, cb_triangles, cb_label, cb_ylabel,
                      cb_tick_interval, cb_minor_tick_interval, cb_offset,
                      cb_width, grid, titlesize, axes_labelsize,
                      tick_labelsize, ax, show, fname])

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
        cb_label : str, optional, default = '$\\theta$ component'
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
        return self.theta.plot(projection=projection,
                               tick_interval=tick_interval,
                               minor_tick_interval=minor_tick_interval,
                               xlabel=xlabel, ylabel=ylabel, title=title,
                               titlesize=titlesize, colorbar=colorbar,
                               cmap=cmap, cmap_limits=cmap_limits,
                               cmap_reverse=cmap_reverse, cb_offset=cb_offset,
                               cb_triangles=cb_triangles, cb_label=cb_label,
                               cb_tick_interval=cb_tick_interval, grid=grid,
                               axes_labelsize=axes_labelsize,
                               cb_ylabel=cb_ylabel, ticks=ticks,
                               cb_minor_tick_interval=cb_minor_tick_interval,
                               tick_labelsize=tick_labelsize, ax=ax,
                               cb_width=cb_width,
                               show=show, fname=fname)

    def plot_phi(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label='$\\phi$ component',
                 cb_tick_interval=None, grid=False, axes_labelsize=None,
                 tick_labelsize=None, show=True, ax=None, cb_offset=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_width=None, fname=None):
        """
        Plot the phi component of the horizontal gradient.

        Usage
        -----
        x.plot_phi([projection, tick_interval, minor_tick_interval, ticks,
                    xlabel, ylabel, title, colorbar, cmap, cmap_limits,
                    cmap_reverse, cb_triangles, cb_label, cb_ylabel,
                    cb_tick_interval, cb_minor_tick_interval, cb_offset,
                    cb_width, grid, titlesize, axes_labelsize, tick_labelsize,
                    ax, show, fname])

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
        cb_label : str, optional, default = '$\\phi$ component'
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
        return self.phi.plot(projection=projection,
                             tick_interval=tick_interval,
                             minor_tick_interval=minor_tick_interval,
                             xlabel=xlabel, ylabel=ylabel, title=title,
                             titlesize=titlesize, colorbar=colorbar,
                             cmap=cmap, cmap_limits=cmap_limits,
                             cmap_reverse=cmap_reverse, cb_offset=cb_offset,
                             cb_triangles=cb_triangles, cb_label=cb_label,
                             cb_tick_interval=cb_tick_interval, grid=grid,
                             axes_labelsize=axes_labelsize,
                             cb_ylabel=cb_ylabel, ticks=ticks,
                             cb_width=cb_width,
                             cb_minor_tick_interval=cb_minor_tick_interval,
                             tick_labelsize=tick_labelsize, ax=ax,
                             show=show, fname=fname)

    def plot(self, projection=None, tick_interval=[60, 30],
             minor_tick_interval=[None, None], xlabel='Longitude',
             ylabel='Latitude', colorbar='bottom', cmap='viridis',
             cmap_limits=None, cmap_reverse=False, cb_triangles='neither',
             cb_tick_interval=None, grid=False, axes_labelsize=9,
             tick_labelsize=8, show=True, cb_offset=None,
             cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
             cb_width=None, fname=None):
        """
        Plot the two vector components of the horizontal gradient.

        Usage
        -----
        x.plot([projection, tick_interval, minor_tick_interval, ticks, xlabel,
                ylabel, colorbar, cmap, cmap_limits, cmap_reverse,
                cb_triangles, cb_ylabel, cb_tick_interval,
                cb_minor_tick_interval, cb_offset, cb_width, grid,
                axes_labelsize, tick_labelsize, show, fname])

        Parameters
        ----------
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        tick_interval : list or tuple, optional, default = [60, 30]
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
        colorbar : str, optional, default = 'bottom'
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
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
        """
        if colorbar is not None:
            if colorbar in set(['bottom', 'top']):
                scale = 0.4
            else:
                scale = 0.25
        else:
            scale = 0.3
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0] * scale)

        fig, ax = _plt.subplots(1, 2, figsize=figsize)
        self.plot_theta(projection=projection, ax=ax.flat[0],
                        tick_interval=tick_interval,
                        minor_tick_interval=minor_tick_interval,
                        xlabel=xlabel, ylabel=ylabel, title=None,
                        titlesize=None, colorbar=colorbar,
                        cmap=cmap, cmap_limits=cmap_limits,
                        cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                        cb_tick_interval=cb_tick_interval,
                        grid=grid, axes_labelsize=axes_labelsize,
                        tick_labelsize=tick_labelsize, cb_offset=cb_offset,
                        cb_ylabel=cb_ylabel, ticks=ticks,
                        cb_minor_tick_interval=cb_minor_tick_interval,
                        cb_width=cb_width, show=show, fname=None)
        self.plot_phi(projection=projection, ax=ax.flat[1],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, title=None,
                      titlesize=None, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      tick_labelsize=tick_labelsize, cb_offset=cb_offset,
                      cb_ylabel=cb_ylabel, ticks=ticks,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, show=show, fname=None)
        fig.tight_layout(pad=0.5)

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def to_xarray(self, title='', description='',
                  comment='pyshtools grid'):
        """
        Return the horizontal gradient gridded data as an xarray DataSet.

        Usage
        -----
        x.to_xarray([title, description, comment])

        Parameters
        ----------
        title : str, optional, default = ''
            Title of the dataset.
        description : str, optional, default = ''
            Description of the dataset ('Remark' in gmt grd files).
        comment : str, optional, default = 'pyshtools grid'
            Additional information about how the data were generated.
        """
        attrs = {'title': title,
                 'description': description,
                 'comment': comment,
                 'nlat': self.nlat,
                 'nlon': self.nlon,
                 'lmax': self.lmax,
                 'grid': self.grid,
                 'lmax_calc': self.lmax_calc,
                 'sampling': self.sampling,
                 'n': self.n,
                 'extend': repr(self.extend)
                 }
        if self.epoch is not None:
            attrs['epoch'] = self.epoch

        _theta = self.theta.to_xarray(title='gradient (theta)',
                                      long_name='theta component',
                                      units=repr(self.units))
        _phi = self.phi.to_xarray(title='gradient (phi)',
                                  long_name='phi components',
                                  units=repr(self.units))

        return _xr.Dataset({'theta': _theta, 'phi': _phi}, attrs=attrs)
