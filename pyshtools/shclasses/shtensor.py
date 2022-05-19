"""
    Class for the gravity and magnetic field 'gradient' tensors.
"""
import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy
from scipy.linalg import eigvalsh as _eigvalsh
import xarray as _xr

from .shgrid import SHGrid as _SHGrid


class Tensor(object):
    """
    Generic class for gravity and magnetic field tensors. To initialize the
    class, use the method tensor() of an SHGravCoeffs or SHMagCoeffs
    class instance.
    """

    def __init__(self):
        """Unused constructor of the main class."""
        print('Initialize the class using one of the two methods:\n'
              '>>> pyshtools.SHGravCoeffs.tensor\n'
              '>>> pyshtools.SHMagCoeffs.tensor\n')

    def compute_invar(self):
        """
        Compute the three invariants (I0, I1, I2) of the tensor, as well as
        the quantity I = -(I2/2)**2 / (I1/3)**3.
        """
        self.i0 = self.vxx + self.vyy + self.vzz
        self.i1 = (self.vxx*self.vyy + self.vyy*self.vzz + self.vxx*self.vzz -
                   self.vxy**2 - self.vyz**2 - self.vxz**2)
        self.i2 = (self.vxx*(self.vyy*self.vzz - self.vyz**2) +
                   self.vxy*(self.vyz*self.vxz - self.vxy*self.vzz) +
                   self.vxz*(self.vxy*self.vyz - self.vxz*self.vyy))
        self.i = (-1.) * (self.i2 / 2.)**2
        self.i.data[1:self.nlat-self.extend, :] /= \
            (self.i1.data[1:self.nlat-self.extend, :] / 3.)**3

    def compute_eig(self):
        """
        Compute the three eigenvalues of the tensor: eig1, eig2, ei3.
        """
        self.eig1 = _SHGrid.from_array(_np.zeros_like(self.vxx.data),
                                       grid='DH')
        self.eig2 = _SHGrid.from_array(_np.zeros_like(self.vxx.data),
                                       grid='DH')
        self.eig3 = _SHGrid.from_array(_np.zeros_like(self.vxx.data),
                                       grid='DH')

        for i in range(self.nlat):
            for j in range(self.nlon):
                a = _np.array([[self.vxx.data[i, j],
                                self.vxy.data[i, j],
                                self.vxz.data[i, j]],
                               [self.vyx.data[i, j],
                                self.vyy.data[i, j],
                                self.vyz.data[i, j]],
                               [self.vzx.data[i, j],
                                self.vzy.data[i, j],
                                self.vzz.data[i, j]]])

                eigs = _eigvalsh(a)

                self.eig1.data[i, j] = eigs[2]
                self.eig2.data[i, j] = eigs[1]
                self.eig3.data[i, j] = eigs[0]

    def compute_eigh(self):
        """
        Compute the two horizontal eigenvalues of the tensor (eigh1, and
        eigh2), as well as the combined maximum absolute value of the two
        (eighh).
        """
        self.eigh1 = _SHGrid.from_array(_np.zeros_like(self.vxx.data),
                                        grid='DH')
        self.eigh2 = _SHGrid.from_array(_np.zeros_like(self.vxx.data),
                                        grid='DH')
        self.eighh = _SHGrid.from_array(_np.zeros_like(self.vxx.data),
                                        grid='DH')

        for i in range(self.nlat):
            for j in range(self.nlon):
                a = _np.array([[self.vxx.data[i, j],
                                self.vxy.data[i, j]],
                               [self.vyx.data[i, j],
                                self.vyy.data[i, j]]])

                eigs = _eigvalsh(a)

                self.eigh1.data[i, j] = eigs[1]
                self.eigh2.data[i, j] = eigs[0]

                if abs(eigs[0]) >= abs(eigs[1]):
                    self.eighh.data[i, j] = eigs[0]
                else:
                    self.eighh.data[i, j] = eigs[1]

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
        Print a summary of the data stored in the SHGravTensor class instance.

        Usage
        -----
        x.info()
        """
        print(repr(self))

    def plot_vxx(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vxx component of the tensor.

        Usage
        -----
        x.plot_vxx([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{xx}$'
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
        if cb_label is None:
            cb_label = self._vxx_label

        return self.vxx.plot(projection=projection,
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

    def plot_vyy(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vyy component of the tensor.

        Usage
        -----
        x.plot_vyy([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{yy}$'
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
        if cb_label is None:
            cb_label = self._vyy_label

        return self.vyy.plot(projection=projection,
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

    def plot_vzz(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vzz component of the tensor.

        Usage
        -----
        x.plot_vzz([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{zz}$'
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
        if cb_label is None:
            cb_label = self._vzz_label

        return self.vzz.plot(projection=projection,
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

    def plot_vxy(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vxx component of the tensor.

        Usage
        -----
        x.plot_vxy([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{xy}$'
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
        if cb_label is None:
            cb_label = self._vxy_label

        return self.vxy.plot(projection=projection,
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

    def plot_vyx(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vyx component of the tensor.

        Usage
        -----
        x.plot_vyx([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{yx}$'
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
        if cb_label is None:
            cb_label = self._vyx_label

        return self.vyx.plot(projection=projection,
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

    def plot_vxz(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vxz component of the tensor.

        Usage
        -----
        x.plot_vxz([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{xz}$'
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
        if cb_label is None:
            cb_label = self._vxz_label

        return self.vxz.plot(projection=projection,
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

    def plot_vzx(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vzx component of the tensor.

        Usage
        -----
        x.plot_vzx([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{zx}$'
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
        if cb_label is None:
            cb_label = self._vzx_label

        return self.vzx.plot(projection=projection,
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

    def plot_vyz(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vyz component of the tensor.

        Usage
        -----
        x.plot_vyz([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{yz}$'
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
        if cb_label is None:
            cb_label = self._vyz_label

        return self.vyz.plot(projection=projection,
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

    def plot_vzy(self, projection=None, tick_interval=[30, 30],
                 minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                 title=None, titlesize=None, colorbar='right',
                 cmap='viridis', cmap_limits=None, cmap_reverse=False,
                 cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                 grid=False, axes_labelsize=None, tick_labelsize=None,
                 cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                 cb_offset=None, cb_width=None, show=True, ax=None,
                 fname=None):
        """
        Plot the Vzy component of the tensor.

        Usage
        -----
        x.plot_vzy([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$V_{zy}$'
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
        if cb_label is None:
            cb_label = self._vzy_label

        return self.vzy.plot(projection=projection,
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

    def plot(self, projection=None, tick_interval=[90, 90],
             minor_tick_interval=[30, 30], xlabel='', ylabel='',
             colorbar='bottom', cmap='viridis', cmap_limits=None,
             cmap_reverse=False, cb_triangles='neither', cb_label=None,
             cb_tick_interval=None, grid=False, axes_labelsize=8,
             cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
             cb_offset=None, cb_width=None, tick_labelsize=8, show=True,
             ax=None, fname=None):
        """
        Plot the 9 components of the tensor.

        Usage
        -----
        x.plot([projection, tick_interval, minor_tick_interval, ticks, xlabel,
                ylabel, colorbar, cmap, cmap_limits, cmap_reverse,
                cb_triangles, cb_label, cb_ylabel, cb_tick_interval,
                cb_minor_tick_interval, cb_offset, cb_width, grid,
                axes_labelsize, tick_labelsize, ax, show, fname])

        Parameters
        ----------
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        tick_interval : list or tuple, optional, default = [90, 90]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
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
        axes_labelsize : int, optional, default = 8
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = 8
            The font size for the x and y tick labels.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
       """
        if colorbar is not None:
            if colorbar in set(['bottom', 'top']):
                scale = 0.9
            else:
                scale = 0.45
        else:
            scale = 0.55
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0] * scale)

        fig, ax = _plt.subplots(3, 3, figsize=figsize)
        self.plot_vxx(projection=projection, ax=ax.flat[0],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vxy(projection=projection, ax=ax.flat[1],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vxz(projection=projection, ax=ax.flat[2],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vyx(projection=projection, ax=ax.flat[3],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vyy(projection=projection, ax=ax.flat[4],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vyz(projection=projection, ax=ax.flat[5],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vzx(projection=projection, ax=ax.flat[6],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vzy(projection=projection, ax=ax.flat[7],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)
        self.plot_vzz(projection=projection, ax=ax.flat[8],
                      tick_interval=tick_interval,
                      minor_tick_interval=minor_tick_interval,
                      xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                      cmap=cmap, cmap_limits=cmap_limits,
                      cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                      cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                      grid=grid, axes_labelsize=axes_labelsize,
                      cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                      cb_minor_tick_interval=cb_minor_tick_interval,
                      cb_width=cb_width, tick_labelsize=tick_labelsize,
                      show=show)

        fig.tight_layout(pad=0.5)

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def plot_i0(self, projection=None, tick_interval=[30, 30],
                minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                title=None, titlesize=None, colorbar='right',
                cmap='viridis', cmap_limits=None, cmap_reverse=False,
                cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                grid=False, axes_labelsize=None, tick_labelsize=None,
                cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                cb_offset=None, cb_width=None, show=True, ax=None, fname=None):
        """
        Plot the first invariant I0 (the trace) of the tensor

            I0 = vxx + vyy + vzz

        which should be identically zero.

        Usage
        -----
        x.plot_i0([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = 'Tr $V_{ij}$'
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
        if cb_label is None:
            cb_label = self._i0_label

        if self.i0 is None:
            self.compute_invar()

        return self.i0.plot(projection=projection,
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

    def plot_i1(self, projection=None, tick_interval=[30, 30],
                minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                title=None, titlesize=None, colorbar='right',
                cmap='viridis', cmap_limits=None, cmap_reverse=False,
                cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                grid=False, axes_labelsize=None, tick_labelsize=None,
                cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                cb_offset=None, cb_width=None, show=True, ax=None, fname=None):
        """
        Plot the second invariant I1 of the tensor:

            I1 = vxx*vyy + vyy*vzz + vxx*vzz - vxy**2 - vyz**2 - vxz**2

        Usage
        -----
        x.plot_i1([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = '$I_1$'
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
        if cb_label is None:
            cb_label = self._i1_label

        if self.i1 is None:
            self.compute_invar()

        return self.i1.plot(projection=projection,
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

    def plot_i2(self, projection=None, tick_interval=[30, 30],
                minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                title=None, titlesize=None, colorbar='right',
                cmap='viridis', cmap_limits=None, cmap_reverse=False,
                cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                grid=False, axes_labelsize=None, tick_labelsize=None,
                cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                cb_offset=None, cb_width=None, show=True, ax=None, fname=None):
        """
        Plot the third invariant I2 (the determinant) of the tensor:

           I2 = vxx*(vyy*vzz - vyz**2) + vxy*(vyz*vxz - vxy*vzz)
                + vxz*(vxy*vyz - vxz*vyy)

        Usage
        -----
        x.plot_i2([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = 'det $V_{ij}$'
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
        if cb_label is None:
            cb_label = self._i2_label

        if self.i2 is None:
            self.compute_invar()

        return self.i2.plot(projection=projection,
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

    def plot_i(self, projection=None, tick_interval=[30, 30],
               minor_tick_interval=[None, None], xlabel=None, ylabel=None,
               title=None, titlesize=None, colorbar='right',
               cmap='viridis', cmap_limits=None, cmap_reverse=False,
               cb_triangles='neither', cb_label=None, cb_tick_interval=None,
               grid=False, axes_labelsize=None, tick_labelsize=None,
               cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
               cb_offset=None, cb_width=None, show=True, ax=None, fname=None):
        """
        Plot the dimensionless quantity I of Pedersen and Rasmussen (1990)

           I = -(I2/2)**2 / (I1/3)**3

        that is bounded by 0 and 1.

        Usage
        -----
        x.plot_i([projection, tick_interval, minor_tick_interval, ticks,
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
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot.
        minor_tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
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
        cb_label : str, optional, default = '$-(I_2/2)^{2} / (I_1/3)^{3}$'
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
        if cb_label is None:
            cb_label = self._i_label

        if self.i is None:
            self.compute_invar()

        return self.i.plot(projection=projection,
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

    def plot_invar(self, projection=None, tick_interval=[60, 60],
                   minor_tick_interval=[30, 30], xlabel='',
                   ylabel='', colorbar='bottom', cmap='viridis',
                   cmap_limits=None, cmap_reverse=False,
                   cb_triangles='neither', cb_label=None,
                   cb_tick_interval=None, grid=False, axes_labelsize=9,
                   cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                   cb_offset=None, cb_width=None, tick_labelsize=8, show=True,
                   ax=None, fname=None):
        """
        Plot the three invariants of the tensor and the derived quantity I.

        Usage
        -----
        x.plot_invar([projection, tick_interval, minor_tick_interval, ticks,
                      xlabel, ylabel, colorbar, cmap, cmap_limits,
                      cmap_reverse, cb_triangles, cb_label, cb_ylabel,
                      cb_tick_interval, cb_minor_tick_interval, cb_offset,
                      cb_width, grid, axes_labelsize, tick_labelsize, ax, show,
                      fname])

        Parameters
        ----------
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        tick_interval : list or tuple, optional, default = [60, 60]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
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
        axes_labelsize : int, optional, default = 9
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = 8
            The font size for the x and y tick labels.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
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

        self.plot_i0(projection=projection, ax=ax.flat[0],
                     tick_interval=tick_interval,
                     minor_tick_interval=minor_tick_interval,
                     xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                     cmap=cmap, cmap_limits=cmap_limits,
                     cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                     cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                     grid=grid, axes_labelsize=axes_labelsize,
                     cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                     cb_minor_tick_interval=cb_minor_tick_interval,
                     cb_width=cb_width, tick_labelsize=tick_labelsize,
                     show=show)
        self.plot_i1(projection=projection, ax=ax.flat[1],
                     tick_interval=tick_interval, cb_offset=cb_offset,
                     minor_tick_interval=minor_tick_interval,
                     xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                     cmap=cmap, cmap_limits=cmap_limits,
                     cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                     cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                     grid=grid, axes_labelsize=axes_labelsize,
                     cb_ylabel=cb_ylabel, ticks=ticks,
                     cb_minor_tick_interval=cb_minor_tick_interval,
                     cb_width=cb_width, tick_labelsize=tick_labelsize,
                     show=show)
        self.plot_i2(projection=projection, ax=ax.flat[2],
                     tick_interval=tick_interval, cb_offset=cb_offset,
                     minor_tick_interval=minor_tick_interval,
                     xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                     cmap=cmap, cmap_limits=cmap_limits,
                     cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                     cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                     grid=grid, axes_labelsize=axes_labelsize,
                     cb_ylabel=cb_ylabel, ticks=ticks,
                     cb_minor_tick_interval=cb_minor_tick_interval,
                     cb_width=cb_width, tick_labelsize=tick_labelsize,
                     show=show)
        self.plot_i(projection=projection, ax=ax.flat[3],
                    tick_interval=tick_interval, cb_offset=cb_offset,
                    minor_tick_interval=minor_tick_interval,
                    xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                    cmap=cmap, cmap_limits=cmap_limits,
                    cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                    cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                    grid=grid, axes_labelsize=axes_labelsize,
                    cb_ylabel=cb_ylabel, ticks=ticks,
                    cb_minor_tick_interval=cb_minor_tick_interval,
                    cb_width=cb_width, tick_labelsize=tick_labelsize,
                    show=show)

        fig.tight_layout(pad=0.5)

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def plot_eig1(self, projection=None, tick_interval=[30, 30],
                  minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                  title=None, titlesize=None, colorbar='right',
                  cmap='viridis', cmap_limits=None, cmap_reverse=False,
                  cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                  grid=False, axes_labelsize=None, tick_labelsize=None,
                  cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                  cb_offset=None, cb_width=None, show=True, ax=None,
                  fname=None):
        r"""
        Plot the first eigenvalue of the tensor.

        Usage
        -----
        x.plot_eig1([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = r'$\lambda_1$'
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
        if cb_label is None:
            cb_label = self._eig1_label

        if self.eig1 is None:
            self.compute_eig()

        return self.eig1.plot(projection=projection,
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

    def plot_eig2(self, projection=None, tick_interval=[30, 30],
                  minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                  title=None, titlesize=None, colorbar='right',
                  cmap='viridis', cmap_limits=None, cmap_reverse=False,
                  cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                  grid=False, axes_labelsize=None, tick_labelsize=None,
                  cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                  cb_offset=None, cb_width=None, show=True, ax=None,
                  fname=None):
        r"""
        Plot the second eigenvalue of the tensor.

        Usage
        -----
        x.plot_eig2([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = r'$\lambda_2$'
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
        if cb_label is None:
            cb_label = self._eig2_label

        if self.eig1 is None:
            self.compute_eig()

        return self.eig2.plot(projection=projection,
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

    def plot_eig3(self, projection=None, tick_interval=[30, 30],
                  minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                  title=None, titlesize=None, colorbar='right',
                  cmap='viridis', cmap_limits=None, cmap_reverse=False,
                  cb_triangles='neither', cb_label=None, cb_tick_interval=None,
                  grid=False, axes_labelsize=None, tick_labelsize=None,
                  cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                  cb_offset=None, cb_width=None, show=True, ax=None,
                  fname=None):
        r"""
        Plot the third eigenvalue of the tensor.

        Usage
        -----
        x.plot_eig3([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = r'$\lambda_3$'
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
        if cb_label is None:
            cb_label = self._eig3_label

        if self.eig1 is None:
            self.compute_eig()

        return self.eig3.plot(projection=projection,
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

    def plot_eigs(self, projection=None, tick_interval=[60, 60],
                  minor_tick_interval=[30, 30], xlabel='',
                  ylabel='', colorbar='bottom', cmap='viridis',
                  cmap_limits=None, cmap_reverse=False,
                  cb_triangles='neither', cb_label=None,
                  cb_tick_interval=None, grid=False, axes_labelsize=9,
                  cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                  cb_offset=None, cb_width=None, tick_labelsize=8, show=True,
                  ax=None, fname=None):
        """
        Plot the three eigenvalues of the tensor.

        Usage
        -----
        x.plot_eigs([projection, tick_interval, minor_tick_interval, ticks,
                     xlabel, ylabel, colorbar, cmap, cmap_limits, cmap_reverse,
                     cb_triangles, cb_label, cb_ylabel, cb_tick_interval,
                     cb_minor_tick_interval, cb_offset, cb_width, grid,
                     axes_labelsize, tick_labelsize, ax, show, fname])

        Parameters
        ----------
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        tick_interval : list or tuple, optional, default = [60, 60]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
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
        axes_labelsize : int, optional, default = 9
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = 8
            The font size for the x and y tick labels.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
        """
        if colorbar is not None:
            if colorbar in set(['bottom', 'top']):
                scale = 2.3
            else:
                scale = 1.4
        else:
            scale = 1.65
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0] * scale)

        fig, ax = _plt.subplots(3, 1, figsize=figsize)

        self.plot_eig1(projection=projection, ax=ax.flat[0],
                       tick_interval=tick_interval,
                       minor_tick_interval=minor_tick_interval,
                       xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                       cmap=cmap, cmap_limits=cmap_limits,
                       cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                       cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                       grid=grid, axes_labelsize=axes_labelsize,
                       cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                       cb_minor_tick_interval=cb_minor_tick_interval,
                       cb_width=cb_width, tick_labelsize=tick_labelsize,
                       show=show)
        self.plot_eig2(projection=projection, ax=ax.flat[1],
                       tick_interval=tick_interval,
                       minor_tick_interval=minor_tick_interval,
                       xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                       cmap=cmap, cmap_limits=cmap_limits,
                       cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                       cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                       grid=grid, axes_labelsize=axes_labelsize,
                       cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                       cb_minor_tick_interval=cb_minor_tick_interval,
                       cb_width=cb_width, tick_labelsize=tick_labelsize,
                       show=show)
        self.plot_eig3(projection=projection, ax=ax.flat[2],
                       tick_interval=tick_interval,
                       minor_tick_interval=minor_tick_interval,
                       xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                       cmap=cmap, cmap_limits=cmap_limits,
                       cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                       cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                       grid=grid, axes_labelsize=axes_labelsize,
                       cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                       cb_minor_tick_interval=cb_minor_tick_interval,
                       cb_width=cb_width, tick_labelsize=tick_labelsize,
                       show=show)

        fig.tight_layout(pad=0.5)

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def plot_eigh1(self, projection=None, tick_interval=[30, 30],
                   minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                   title=None, titlesize=None, colorbar='right',
                   cmap='viridis', cmap_limits=None, cmap_reverse=False,
                   cb_triangles='neither', cb_label=None,
                   cb_tick_interval=None, grid=False, axes_labelsize=None,
                   cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                   cb_offset=None, cb_width=None, tick_labelsize=None,
                   show=True, ax=None, fname=None):
        r"""
        Plot the first eigenvalue of the horizontal tensor.

        Usage
        -----
        x.plot_eigh1([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = r'$\lambda_{h1}$'
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
        if cb_label is None:
            cb_label = self._eigh1_label

        if self.eigh1 is None:
            self.compute_eigh()

        return self.eigh1.plot(projection=projection,
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

    def plot_eigh2(self, projection=None, tick_interval=[30, 30],
                   minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                   title=None, titlesize=None, colorbar='right',
                   cmap='viridis', cmap_limits=None, cmap_reverse=False,
                   cb_triangles='neither', cb_label=None,
                   cb_tick_interval=None, grid=False, axes_labelsize=None,
                   cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                   cb_offset=None, cb_width=None, tick_labelsize=None,
                   show=True, ax=None, fname=None):
        r"""
        Plot the second eigenvalue of the horizontal tensor.

        Usage
        -----
        x.plot_eigh2([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = r'$\lambda_{h2}$'
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
        if cb_label is None:
            cb_label = self._eigh2_label

        if self.eigh1 is None:
            self.compute_eigh()

        return self.eigh2.plot(projection=projection,
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

    def plot_eighh(self, projection=None, tick_interval=[30, 30],
                   minor_tick_interval=[None, None], xlabel=None, ylabel=None,
                   title=None, titlesize=None, colorbar='right',
                   cmap='viridis', cmap_limits=None, cmap_reverse=False,
                   cb_triangles='neither', cb_label=None,
                   cb_tick_interval=None, grid=False, axes_labelsize=None,
                   cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                   cb_offset=None, cb_width=None, tick_labelsize=None,
                   show=True, ax=None, fname=None):
        r"""
        Plot the maximum absolute value eigenvalue of the horizontal tensor.

        Usage
        -----
        x.plot_eighh([projection, tick_interval, minor_tick_interval, ticks,
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
        cb_label : str, optional, default = r'$\lambda_{hh}$'
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
        if cb_label is None:
            cb_label = self._eighh_label

        if self.eigh1 is None:
            self.compute_eigh()

        return self.eighh.plot(projection=projection,
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

    def plot_eigh(self, projection=None, tick_interval=[60, 60],
                  minor_tick_interval=[30, 30], xlabel='',
                  ylabel='', colorbar='bottom', cmap='viridis',
                  cmap_limits=None, cmap_reverse=False,
                  cb_triangles='neither', cb_label=None,
                  cb_tick_interval=None, grid=False, axes_labelsize=9,
                  cb_minor_tick_interval=None, ticks='WSen', cb_ylabel=None,
                  cb_offset=None, cb_width=None, tick_labelsize=8, show=True,
                  ax=None, fname=None):
        """
        Plot the two eigenvalues and maximum absolute value eigenvalue of the
        horizontal tensor.

        Usage
        -----
        x.plot_eigh([projection, tick_interval, minor_tick_interval, ticks,
                     xlabel, ylabel, colorbar, cmap, cmap_limits, cmap_reverse,
                     cb_triangles, cb_label, cb_ylabel, cb_tick_interval,
                     cb_minor_tick_interval, cb_offset, cb_width, grid,
                     axes_labelsize, tick_labelsize, ax, show, fname])

        Parameters
        ----------
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        tick_interval : list or tuple, optional, default = [60, 60]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
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
        axes_labelsize : int, optional, default = 9
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = 8
            The font size for the x and y tick labels.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
        """
        if colorbar is not None:
            if colorbar in set(['bottom', 'top']):
                scale = 2.3
            else:
                scale = 1.4
        else:
            scale = 1.65
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0] * scale)

        fig, ax = _plt.subplots(3, 1, figsize=figsize)

        self.plot_eigh1(projection=projection, ax=ax.flat[0],
                        tick_interval=tick_interval,
                        minor_tick_interval=minor_tick_interval,
                        xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                        cmap=cmap, cmap_limits=cmap_limits,
                        cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                        cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                        grid=grid, axes_labelsize=axes_labelsize,
                        cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                        cb_minor_tick_interval=cb_minor_tick_interval,
                        cb_width=cb_width, tick_labelsize=tick_labelsize,
                        show=show)
        self.plot_eigh2(projection=projection, ax=ax.flat[1],
                        tick_interval=tick_interval,
                        minor_tick_interval=minor_tick_interval,
                        xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                        cmap=cmap, cmap_limits=cmap_limits,
                        cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                        cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                        grid=grid, axes_labelsize=axes_labelsize,
                        cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                        cb_minor_tick_interval=cb_minor_tick_interval,
                        cb_width=cb_width, tick_labelsize=tick_labelsize,
                        show=show)
        self.plot_eighh(projection=projection, ax=ax.flat[2],
                        tick_interval=tick_interval,
                        minor_tick_interval=minor_tick_interval,
                        xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                        cmap=cmap, cmap_limits=cmap_limits,
                        cmap_reverse=cmap_reverse, cb_triangles=cb_triangles,
                        cb_label=cb_label, cb_tick_interval=cb_tick_interval,
                        grid=grid, axes_labelsize=axes_labelsize,
                        cb_ylabel=cb_ylabel, ticks=ticks, cb_offset=cb_offset,
                        cb_minor_tick_interval=cb_minor_tick_interval,
                        cb_width=cb_width, tick_labelsize=tick_labelsize,
                        show=show)

        fig.tight_layout(pad=0.5)

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def to_xarray(self, title='', description='',
                  comment='pyshtools grid'):
        """
        Return all tensor gridded data as an xarray DataSet.

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
                 'lmax_calc': self.lmax_calc,
                 'sampling': self.sampling,
                 'grid': self.grid,
                 'a': self.a,
                 'f': self.f,
                 'n': self.n,
                 'extend': repr(self.extend)
                 }
        if isinstance(self, SHGravTensor):
            attrs['gm'] = self.gm
            if self.epoch is not None:
                attrs['epoch'] = self.epoch
            desc = 'gravity tensor component '
        else:
            if self.year is not None:
                attrs['year'] = self.year
            desc = 'magnetic field tensor component '

        _vxx = self.vxx.to_xarray(title=desc+'(Vxx)', long_name='$V_{xx}$',
                                  units=self._vii_units)
        _vxy = self.vxy.to_xarray(title=desc+'(Vxy)', long_name='$V_{xy}$',
                                  units=self._vii_units)
        _vxz = self.vxz.to_xarray(title=desc+'(Vxz)', long_name='$V_{xz}$',
                                  units=self._vii_units)
        _vyx = self.vyx.to_xarray(title=desc+'(Vyx)', long_name='$V_{yx}$',
                                  units=self._vii_units)
        _vyy = self.vyy.to_xarray(title=desc+'(Vyy)', long_name='$V_{yy}$',
                                  units=self._vii_units)
        _vyz = self.vyz.to_xarray(title=desc+'(Vyz)', long_name='$V_{yz}$',
                                  units=self._vii_units)
        _vzx = self.vzx.to_xarray(title=desc+'(Vzx)', long_name='$V_{zx}$',
                                  units=self._vii_units)
        _vzy = self.vzy.to_xarray(title=desc+'(Vzy)', long_name='$V_{zy}$',
                                  units=self._vii_units)
        _vzz = self.vzz.to_xarray(title=desc+'(Vzz)', long_name='$V_{zz}$',
                                  units=self._vii_units)

        dataset = _xr.Dataset({'vxx': _vxx, 'vxy': _vxy, 'vxz': _vxz,
                               'vyx': _vyx, 'vyy': _vyy, 'vyz': _vyz,
                               'vzx': _vzx, 'vzy': _vzy, 'vzz': _vzz},
                              attrs=attrs)

        if self.i0 is not None:
            if isinstance(self, SHGravTensor):
                desc0 = 'First invariant of the gravity tensor'
                desc1 = 'Second invariant of the gravity tensor'
                desc2 = 'Third invariant of the gravity tensor'
                desc = 'Unitless invariant of the gravity tensor'
            else:
                desc0 = 'First invariant of the magnetic field tensor'
                desc1 = 'Second invariant of the magnetic field tensor'
                desc2 = 'Third invariant of the magnetic field tensor'
                desc = 'Unitless invariant of the magnetic field tensor'

            _i0 = self.i0.to_xarray(title=desc0,
                                    long_name='$I_0$, Tr $V_{ii}$',
                                    units=self._i0_units)
            _i1 = self.i1.to_xarray(title=desc1, long_name='$I_1$',
                                    units=self._i1_units)
            _i2 = self.i2.to_xarray(title=desc2,
                                    long_name='$I_2$, det $V_{ij}$',
                                    units=self._i2_units)
            _i = self.i.to_xarray(title=desc,
                                  long_name='$-(I_2/2)^{2} / ' +
                                  '(I_1/3)^{3}$',
                                  units='none')

            dataset['i0'] = _i0
            dataset['i1'] = _i1
            dataset['i2'] = _i2
            dataset['i'] = _i

        if self.eig1 is not None:
            if isinstance(self, SHGravTensor):
                desc1 = 'First eigenvalue of the gravity tensor'
                desc2 = 'Second eigenvalue of the gravity tensor'
                desc3 = 'Third eigenvalue of the gravity tensor'
            else:
                desc1 = 'First eigenvalue of the magnetic field tensor'
                desc2 = 'Second eigenvalue of the magnetic field tensor'
                desc3 = 'Third eigenvalue of the magnetic field tensor'

            _eig1 = self.eig1.to_xarray(title=desc1,
                                        long_name=r'${\lambda}_1$',
                                        units=self._vii_units)
            _eig2 = self.eig2.to_xarray(title=desc2,
                                        long_name=r'${\lambda}_2$',
                                        units=self._vii_units)
            _eig3 = self.eig3.to_xarray(title=desc3,
                                        long_name=r'${\lambda}_3$',
                                        units=self._vii_units)

            dataset['eig1'] = _eig1
            dataset['eig2'] = _eig2
            dataset['eig3'] = _eig3

        if self.eighh is not None:
            if isinstance(self, SHGravTensor):
                desc1 = 'First horizontal eigenvalue of the gravity tensor'
                desc2 = 'Second horizontal eigenvalue of the gravity tensor'
                desc3 = 'Combined horizontal eigenvalue of the gravity tensor'
            else:
                desc1 = 'First horizontal eigenvalue of the magnetic ' \
                        + 'field tensor'
                desc2 = 'Second horizontal eigenvalue of the magnetic ' \
                        + 'field tensor'
                desc3 = 'Combined horizontal eigenvalue of the magnetic ' \
                        + 'field tensor'

            _eigh1 = self.eigh1.to_xarray(title=desc1,
                                          long_name=r'${\lambda}_{h1}$',
                                          units=self._vii_units)
            _eigh2 = self.eigh2.to_xarray(title=desc2,
                                          long_name=r'${\lambda}_{h2}$',
                                          units=self._vii_units)
            _eighh = self.eighh.to_xarray(title=desc3,
                                          long_name=r'${\lambda}_{hh}$',
                                          units=self._vii_units)

            dataset['eigh1'] = _eigh1
            dataset['eigh2'] = _eigh2
            dataset['eighh'] = _eighh

        return dataset


class SHGravTensor(Tensor):
    """
    Class for the gravity field tensor and eigenvalues. The class is
    initialized from a class instance of SHGravCoeffs using the method
    tensor().

    Attributes:

    vxx, vxy, vzz,   : The 9 components of the gravity tensor.
    vyx, vyy, vyz,
    vzx, vzy, vzz
    i0, i1, i2, i    : The three invariants of the gravity tensor and a
                       derived quantity that is bounded between 0 and 1.
                       These are computed by a call to compute_invar().
    eig1, eig2, eig3 : The three eigenvalues of the gravity tensor, which are
                       computed by a call to compute_eig().
    eigh1, eigh2,    : The horizontal eigenvalues of the gravity tensor, which
    eighh              are computed by a call to compute_eigh().
    gm               : The gravitational constant times the mass of the body.
    a                : Semimajor axis of the reference ellipsoid.
    f                : Flattening of the reference ellipsoid, f=(a-b)/a.
    lmax             : The maximum spherical harmonic degree resolvable by the
                       grids.
    lmax_calc        : The maximum spherical harmonic degree of the
                       gravitational potential used in creating the grids.
    units            : The units of the gridded data.
    epoch          : The epoch time of the gravity model.
    nlat, nlon       : The number of latitude and longitude bands in the grids.
    n                : The number of samples in latitude.
    sampling         : The longitudinal sampling for Driscoll and Healy grids.
                       Either 1 for equally sampled grids (nlat=nlon) or 2 for
                       equally spaced grids in degrees.
    extend           : True if the grid contains the redundant column for
                       360 E and the unnecessary row for 90 S.

    Methods:

    plot()          : Plot all 9 components of the gravity tensor.
    plot_vxx()      : Plot the vxx component of the gravity tensor.
    plot_vxy()      : Plot the vxy component of the gravity tensor.
    plot_vxz()      : Plot the vxz component of the gravity tensor.
    plot_vyx()      : Plot the vyx component of the gravity tensor.
    plot_vyy()      : Plot the vyy component of the gravity tensor.
    plot_vyz()      : Plot the vyz component of the gravity tensor.
    plot_vzx()      : Plot the vzx component of the gravity tensor.
    plot_vzy()      : Plot the vzy component of the gravity tensor.
    plot_vzz()      : Plot the vzz component of the gravity tensor.

    compute_invar() : Compute the invariants of the gravity tensor.
    plot_i0()       : Plot the first invariant I0 of the gravity tensor.
    plot_i1()       : Plot the second invariant I1 of the gravity tensor.
    plot_i2()       : Plot the third invariant I2 of the gravity tensor.
    plot_i()        : Plot the derived quantity I = -(I2/2)**2 / (I1/3)**3.

    compute_eig()   : Compute the three eigenvalues of the gravity tensor.
    plot_eig()      : Plot the three eigenvalues of the gravity tensor.
    plot_eig1()     : Plot the first eigenvalue of the gravity tensor.
    plot_eig2()     : Plot the second eigenvalue of the gravity tensor.
    plot_eig3()     : Plot the third eigenvalue of the gravity tensor.

    compute_eigh()  : Compute the horizontal eigenvalues of the gravity tensor.
    plot_eigh()     : Plot the two horizontal eigenvalues and the combined
                      maximum absolute eigenvalue of the gravity tensor.
    plot_eigh1()    : Plot the first horizontal eigenvalue of the gravity
                      tensor.
    plot_eigh2()    : Plot the second horizontal eigenvalue of the gravity
                      tensor.
    plot_eighh()    : Plot the combined maximum absolute eigenvalue of the
                      gravity tensor.

    to_xarray()     : Return an xarray DataSet of all gridded data.
    copy()          : Return a copy of the class instance.
    info()          : Print a summary of the data stored in the SHGravTensor
                      instance.
    """

    def __init__(self, vxx, vyy, vzz, vxy, vxz, vyz, gm, a, f, lmax,
                 lmax_calc, units='Eötvös', epoch=None):
        """
        Initialize the SHGravTensor class.
        """
        self.vxx = _SHGrid.from_array(vxx, grid='DH', units=units)
        self.vyy = _SHGrid.from_array(vyy, grid='DH', units=units)
        self.vzz = _SHGrid.from_array(vzz, grid='DH', units=units)
        self.vxy = _SHGrid.from_array(vxy, grid='DH', units=units)
        self.vxz = _SHGrid.from_array(vxz, grid='DH', units=units)
        self.vyz = _SHGrid.from_array(vyz, grid='DH', units=units)
        self.vyx = self.vxy
        self.vzx = self.vxz
        self.vzy = self.vyz
        self.grid = self.vxx.grid
        self.sampling = self.vxx.sampling
        self.nlat = self.vxx.nlat
        self.nlon = self.vxx.nlon
        self.n = self.vxx.n
        self.extend = self.vxx.extend
        self.gm = gm
        self.a = a
        self.f = f
        self.lmax = lmax
        self.lmax_calc = lmax_calc
        self.i0 = None
        self.i1 = None
        self.i2 = None
        self.i = None
        self.eig1 = None
        self.eig2 = None
        self.eig3 = None
        self.eigh1 = None
        self.eigh2 = None
        self.eighh = None
        self.units = units
        self.epoch = epoch

        self._vxx_label = '$V_{xx}$, ' + self.units
        self._vxy_label = '$V_{xy}$, ' + self.units
        self._vxz_label = '$V_{xz}$, ' + self.units
        self._vyx_label = '$V_{yx}$, ' + self.units
        self._vyy_label = '$V_{yy}$, ' + self.units
        self._vyz_label = '$V_{yz}$, ' + self.units
        self._vzx_label = '$V_{zx}$, ' + self.units
        self._vzy_label = '$V_{zy}$, ' + self.units
        self._vzz_label = '$V_{zz}$, ' + self.units

        self._i0_label = 'Tr $V_{ii}$, ' + self.units
        self._i1_label = '$I_1$, ' + self.units + '$^2$'
        self._i2_label = 'det $V_{ij}$, ' + self.units + '$^3$'
        self._i_label = '$-(I_2/2)^{2} / (I_1/3)^{3}$'

        self._eig1_label = r'$\lambda_1$, ' + self.units
        self._eig2_label = r'$\lambda_2$, ' + self.units
        self._eig3_label = r'$\lambda_3$, ' + self.units

        self._eigh1_label = r'$\lambda_{h1}$, ' + self.units
        self._eigh2_label = r'$\lambda_{h2}$, ' + self.units
        self._eighh_label = r'$\lambda_{hh}$, ' + self.units

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
               'units = {:s}\n'
               'epoch = {:s}'
               .format(self.grid, self.nlat, self.nlon, self.n, self.sampling,
                       self.extend, self.lmax, self.lmax_calc, self.gm, self.a,
                       self.f, repr(self.units), repr(self.epoch)))
        return str


class SHMagTensor(Tensor):
    """
    Class for the magnetic field tensor and eigenvalues. The class is
    initialized from a class instance of SHMagCoeffs using the method
    tensor().

    Attributes:

    vxx, vxy, vzz,   : The 9 components of the magnetic field tensor.
    vyx, vyy, vyz,
    vzx, vzy, vzz
    i0, i1, i2, i    : The three invariants of the magnetic field tensor and a
                       derived quantity that is bounded between 0 and 1.
    eig1, eig2, eig3 : The three eigenvalues of the magnetic field tensor,
                       which are computed by a call to compute_eig().
    eigh1, eigh2,    : The horizontal eigenvalues of the magnetic field
    eighh              tensor, which are computed by a call to compute_eigh().
    a                : Semimajor axis of the reference ellipsoid.
    f                : Flattening of the reference ellipsoid, f=(a-b)/a.
    lmax             : The maximum spherical harmonic degree resolvable by the
                       grids.
    lmax_calc        : The maximum spherical harmonic degree of the
                       magnetic potential used in creating the grids.
    units            : The units of the gridded data.
    year             : The year of the time-variable magnetic field data.
    nlat, nlon       : The number of latitude and longitude bands in the grids.
    sampling         : The longitudinal sampling for Driscoll and Healy grids.
                       Either 1 for equally sampled grids (nlat=nlon) or 2 for
                       equally spaced grids in degrees.
    extend           : True if the grid contains the redundant column for
                       360 E and the unnecessary row for 90 S.

    Methods:

    plot()          : Plot all 9 components of the magnetic field tensor.
    plot_vxx()      : Plot the vxx component of the magnetic field tensor.
    plot_vxy()      : Plot the vxy component of the magnetic field tensor.
    plot_vxz()      : Plot the vxz component of the magnetic field tensor.
    plot_vyx()      : Plot the vyx component of the magnetic field tensor.
    plot_vyy()      : Plot the vyy component of the magnetic field tensor.
    plot_vyz()      : Plot the vyz component of the magnetic field tensor.
    plot_vzx()      : Plot the vzx component of the magnetic field tensor.
    plot_vzy()      : Plot the vzy component of the magnetic field tensor.
    plot_vzz()      : Plot the vzz component of the magnetic field tensor.

    compute_invar() : Compute the invariants of the magnetic field tensor.
    plot_i0()       : Plot the first invariant I0 of the magnetic field tensor.
    plot_i1()       : Plot the second invariant I1 of themagnetic field tensor.
    plot_i2()       : Plot the third invariant I2 of the magnetic field tensor.
    plot_i()        : Plot the derived quantity I = -(I2/2)**2 / (I1/3)**3.

    compute_eig()   : Compute the three eigenvalues of the magnetic field
                      tensor.
    plot_eig()      : Plot the three eigenvalues of the magnetic field tensor.
    plot_eig1()     : Plot the first eigenvalue of the magnetic field tensor.
    plot_eig2()     : Plot the second eigenvalue of the magnetic field tensor.
    plot_eig3()     : Plot the third eigenvalue of the magnetic field tensor.

    compute_eigh()  : Compute the horizontal eigenvalues of the magnetic field
                      tensor.
    plot_eigh()     : Plot the two horizontal eigenvalues and the combined
                      maximum absolute eigenvalue of the magnetic field tensor.
    plot_eigh1()    : Plot the first horizontal eigenvalue of the magnetic
                      field tensor.
    plot_eigh2()    : Plot the second horizontal eigenvalue of the magnetic
                      field tensor.
    plot_eighh()    : Plot the combined maximum absolute eigenvalue of the
                      magnetic field tensor.

    to_xarray()     : Return an xarray DataSet of all gridded data.
    copy()          : Return a copy of the class instance.
    info()          : Print a summary of the data stored in the SHMagTensor
                      instance.
    """

    def __init__(self, vxx, vyy, vzz, vxy, vxz, vyz, a, f, lmax,
                 lmax_calc, units=None, year=None):
        """
        Initialize the SHMagTensor class.
        """
        self.vxx = _SHGrid.from_array(vxx, grid='DH', units=units)
        self.vyy = _SHGrid.from_array(vyy, grid='DH', units=units)
        self.vzz = _SHGrid.from_array(vzz, grid='DH', units=units)
        self.vxy = _SHGrid.from_array(vxy, grid='DH', units=units)
        self.vxz = _SHGrid.from_array(vxz, grid='DH', units=units)
        self.vyz = _SHGrid.from_array(vyz, grid='DH', units=units)
        self.vyx = self.vxy
        self.vzx = self.vxz
        self.vzy = self.vyz
        self.grid = self.vxx.grid
        self.sampling = self.vxx.sampling
        self.nlat = self.vxx.nlat
        self.nlon = self.vxx.nlon
        self.n = self.vxx.n
        self.extend = self.vxx.extend
        self.a = a
        self.f = f
        self.lmax = lmax
        self.lmax_calc = lmax_calc
        self.i0 = None
        self.i1 = None
        self.i2 = None
        self.i = None
        self.eig1 = None
        self.eig2 = None
        self.eig3 = None
        self.eigh1 = None
        self.eigh2 = None
        self.eighh = None
        self.units = units
        self.year = year
        if self.units.lower() == 'nt/m':
            self._units_formatted = 'nT m$^{-1}$'
            self._i1_units = 'nT$^2$ m$^{-2}$'
            self._i2_units = 'nT$^3$ m$^{-3}$'
        else:
            self._units_formatted = 'T m$^{-1}$'
            self._i1_units = 'T$^2$ m$^{-2}$'
            self._i2_units = 'T$^3$ m$^{-3}$'

        self._vxx_label = '$V_{xx}$, ' + self._units_formatted
        self._vxy_label = '$V_{xy}$, ' + self._units_formatted
        self._vxz_label = '$V_{xz}$, ' + self._units_formatted
        self._vyx_label = '$V_{yx}$, ' + self._units_formatted
        self._vyy_label = '$V_{yy}$, ' + self._units_formatted
        self._vyz_label = '$V_{yz}$, ' + self._units_formatted
        self._vzx_label = '$V_{zx}$, ' + self._units_formatted
        self._vzy_label = '$V_{zy}$, ' + self._units_formatted
        self._vzz_label = '$V_{zz}$, ' + self._units_formatted

        self._i0_label = 'Tr $V_{ii}$, ' + self._units_formatted
        self._i1_label = '$I_1$, ' + self._i1_units
        self._i2_label = 'det $V_{ij}$, ' + self._i2_units
        self._i_label = '$-(I_2/2)^{2} / (I_1/3)^{3}$'

        self._eig1_label = r'$\lambda_1$, ' + self._units_formatted
        self._eig2_label = r'$\lambda_2$, ' + self._units_formatted
        self._eig3_label = r'$\lambda_3$, ' + self._units_formatted

        self._eigh1_label = r'$\lambda_{h1}$, ' + self._units_formatted
        self._eigh2_label = r'$\lambda_{h2}$, ' + self._units_formatted
        self._eighh_label = r'$\lambda_{hh}$, ' + self._units_formatted

    def __repr__(self):
        str = ('grid = {:s}\n'
               'nlat = {:d}\n'
               'nlon = {:d}\n'
               'n = {:d}\n'
               'sampling = {:d}\n'
               'extend = {}\n'
               'lmax = {:d}\n'
               'lmax_calc = {:d}\n'
               'a (m)= {:e}\n'
               'f = {:e}\n'
               'units = {:s}\n'
               'year = {:s}'
               .format(self.grid, self.nlat, self.nlon, self.n, self.sampling,
                       self.extend, self.lmax, self.lmax_calc, self.a,
                       self.f, repr(self.units), repr(self.year)))
        return str
