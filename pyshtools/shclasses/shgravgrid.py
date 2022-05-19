"""
    Class for grids of the three components of the gravity field, the
    gravitational disturbance, and the gravitational potential.
"""
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy
import xarray as _xr

from .shgrid import SHGrid as _SHGrid


class SHGravGrid(object):
    """
    Class for grids of the gravitational potential, three vector components of
    the gravity field, and the total gravitational disturbance. The class is
    initialized from a class instance of SHGravCoeffs using the method
    expand().

    Attributes:

    rad            : SHGrid class instance of the radial component of the
                     gravitational acceleration evaluated on an ellipsoid.
    theta          : SHGrid class instance of the theta component of the
                     gravitational acceleration evaluated on an ellipsoid.
    phi            : SHGrid class instance of the phi component of the
                     gravitational acceleration evaluated on an ellipsoid.
    total          : SHGrid class instance of the total gravitational
                     acceleration on an ellipsoid.
    pot            : SHGrid class instance of the gravitational potential
                     evaluated on an ellipsoid.
    gm             : Gravitational constant times the mass of the body.
    a              : Semimajor axis of the reference ellipsoid.
    f              : Flattening of the reference ellipsoid, f=(a-b)/a.
    omega          : Angular rotation rate of the body.
    normal_gravity : True if the normal gravity is removed from the total
                     gravitational acceleration.
    units          : The units of the gridded gravity data.
    pot_units      : The units of the gridded gravitational potential data.
    epoch          : The epoch time of the gravity model.
    lmax           : The maximum spherical harmonic degree resolvable by the
                     grids.
    lmax_calc      : The maximum spherical harmonic degree of the gravitational
                     potential used in creating the grid.
    nlat, nlon     : The number of latitude and longitude bands in the grids.
    n              : The number of samples in latitude.
    sampling       : The longitudinal sampling for Driscoll and Healy grids.
                     Either 1 for equally sampled grids (nlat=nlon) or 2 for
                     equally spaced grids in degrees.
    extend         : True if the grid contains the redundant column for 360 E
                     and the unnecessary row for 90 S.


    Methods:

    plot()        : Plot all three components of the gravity field and the
                    total gravity disturbance.
    plot_rad()    : Plot the radial component of the gravity field.
    plot_theta()  : Plot the theta component of the gravity field.
    plot_phi()    : Plot the phi component of the gravity field.
    plot_total()  : Plot the total gravity disturbance.
    plot_pot()    : Plot the gravitational potential.
    to_xarray()   : Return an xarray DataSet of all gridded data.
    copy()        : Return a copy of the class instance.
    info()        : Print a summary of the data stored in the SHGravGrid
                    instance.
    """

    def __init__(self, rad, theta, phi, total, pot, gm, a, f, omega,
                 normal_gravity, lmax, lmax_calc, units=None, pot_units=None,
                 epoch=None):
        """
        Initialize the SHGravGrid class.
        """
        self.rad = _SHGrid.from_array(rad, grid='DH', units=units)
        self.theta = _SHGrid.from_array(theta, grid='DH', units=units)
        self.phi = _SHGrid.from_array(phi, grid='DH', units=units)
        self.total = _SHGrid.from_array(total, grid='DH', units=units)
        self.pot = _SHGrid.from_array(pot, grid='DH', units=pot_units)
        self.grid = self.rad.grid
        self.sampling = self.rad.sampling
        self.nlat = self.rad.nlat
        self.nlon = self.rad.nlon
        self.n = self.rad.n
        self.extend = self.rad.extend
        self.gm = gm
        self.a = a
        self.f = f
        if omega is None:
            self.omega = 0.0
        else:
            self.omega = omega
        self.normal_gravity = normal_gravity
        self.lmax = lmax
        self.lmax_calc = lmax_calc
        self.units = units
        self.pot_units = pot_units
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
        Print a summary of the data stored in the SHGravGrid class instance.

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
               'gm (m3 / s2) = {:e}\n'
               'a (m)= {:e}\n'
               'f = {:e}\n'
               'omega (rad / s) = {:s}\n'
               'normal gravity is removed = {:s}\n'
               'units (gravity) = {:s}\n'
               'units (potential) = {:s}\n'
               'epoch = {:s}'
               .format(self.grid, self.nlat, self.nlon, self.n, self.sampling,
                       self.extend, self.lmax, self.lmax_calc, self.gm,
                       self.a, self.f, repr(self.omega),
                       repr(self.normal_gravity), repr(self.units),
                       repr(self.pot_units), repr(self.epoch)))
        return str

    def plot_rad(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label='$g_r$, m s$^{-2}$',
                 cb_tick_interval=None, cb_offset=None,
                 cb_width=None, grid=False, axes_labelsize=None,
                 tick_labelsize=None, show=True, ax=None,
                 cb_minor_tick_interval=None, ticks='WSen',
                 cb_ylabel=None, fname=None):
        """
        Plot the radial component of the gravity field.

        Usage
        -----
        x.plot_rad([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$g_r$, m s$^{-2}$'
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
        return self.rad.plot(projection=projection,
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
                             cb_width=cb_width,
                             tick_labelsize=tick_labelsize, ax=ax,
                             show=show, fname=fname)

    def plot_theta(self, projection=None, tick_interval=[30, 30],
                   minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                   title=None, titlesize=None, colorbar='right',
                   cmap='viridis', cmap_limits=None, cmap_reverse=False,
                   cb_triangles='neither', cb_label=r'$g_\theta$, m s$^{-2}$',
                   cb_tick_interval=None, grid=False, axes_labelsize=None,
                   tick_labelsize=None, show=True, ax=None, cb_offset=None,
                   cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                   fname=None, cb_width=None):
        r"""
        Plot the theta component of the gravity field.

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
        cb_label : str, optional, default = r'$g_\theta$, m s$^{-2}$'
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
                 cb_triangles='neither', cb_label=r'$g_\phi$, m s$^{-2}$',
                 cb_tick_interval=None, grid=False, axes_labelsize=None,
                 tick_labelsize=None, show=True, ax=None, cb_offset=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_width=None, fname=None):
        r"""
        Plot the phi component of the gravity field.

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
        cb_label : str, optional, default = r'$g_\phi$, m s$^{-2}$'
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

    def plot_total(self, projection=None, tick_interval=[30, 30],
                   minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                   title=None, titlesize=None, colorbar='right',
                   cmap='viridis', cmap_limits=None, cmap_reverse=False,
                   cb_triangles='neither', cb_label='Gravity disturbance',
                   cb_tick_interval=None, grid=False, axes_labelsize=None,
                   tick_labelsize=None, show=True, ax=None, cb_offset=None,
                   cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                   cb_width=None, fname=None):
        """
        Plot the total gravity disturbance.

        Usage
        -----
        x.plot_total([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = 'Gravity disturbance'
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

        Note
        -----
        If the normal gravity is removed from the total gravitational
        acceleration, the output will be displayed in mGals.
        """
        if self.normal_gravity:
            if cb_label is not None and cb_ylabel is None:
                cb_label += ', mGal'
        else:
            if cb_label is not None and cb_ylabel is None:
                cb_label += ', m s$^{-2}$'

        if self.normal_gravity:
            temp = self.total*1.e5
        else:
            temp = self.total

        return temp.plot(projection=projection,
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

    def plot_pot(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label='Potential, m$^2$ s$^{-2}$',
                 cb_tick_interval=None, grid=False, axes_labelsize=None,
                 tick_labelsize=None, show=True, ax=None, cb_offset=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_width=None, fname=None):
        """
        Plot the gravitational potential.

        Usage
        -----
        x.plot_pot([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = 'Potential, m$^2$ s$^{-2}$'
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
        return self.pot.plot(projection=projection,
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
                             cb_width=cb_width,
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
        Plot the three vector components of the gravity field and the gravity
        disturbance.

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
                scale = 0.8
            else:
                scale = 0.5
        else:
            scale = 0.6
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0] * scale)

        fig, ax = _plt.subplots(2, 2, figsize=figsize)
        self.plot_rad(projection=projection, ax=ax.flat[0],
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
        self.plot_theta(projection=projection, ax=ax.flat[1],
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
        self.plot_phi(projection=projection, ax=ax.flat[2],
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
        self.plot_total(projection=projection, ax=ax.flat[3],
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
        Return the gravity gridded data as an xarray DataSet.

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
                 'gm': self.gm,
                 'a': self.a,
                 'f': self.f,
                 'omega': self.omega,
                 'lmax_calc': self.lmax_calc,
                 'sampling': self.sampling,
                 'normal_gravity': self.normal_gravity,
                 'n': self.n,
                 'extend': repr(self.extend)
                 }
        if self.epoch is not None:
            attrs['epoch'] = self.epoch

        _total = self.total.to_xarray(title='total gravity disturbance',
                                      long_name='$|g|$', units='$m s^{-2}$')
        _rad = self.rad.to_xarray(title='gravity (radial)',
                                  long_name='$g_r$', units='$m s^{-2}$')
        _theta = self.theta.to_xarray(title='gravity (theta)',
                                      long_name=r'$g_\theta$',
                                      units='$m s^{-2}$')
        _phi = self.phi.to_xarray(title='gravity (phi)',
                                  long_name=r'$g_\phi$', units='$m s^{-2}$')
        _pot = self.pot.to_xarray(title='gravitational potential',
                                  long_name='potential', units='$m^2 s^{-2}$')

        return _xr.Dataset({'radial': _rad, 'theta': _theta, 'phi': _phi,
                            'total': _total, 'potential': _pot}, attrs=attrs)
