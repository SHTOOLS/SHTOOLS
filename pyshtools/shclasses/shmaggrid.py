"""
    Class for grids of the three components of the magnetic field, the
    magnetic intensity, and the magnetic potential.
"""
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy

from .shcoeffsgrid import SHGrid as _SHGrid


class SHMagGrid(object):
    """
    Class for grids of the magnetic potential, three vector components of
    the magnetic field, and the total magnetic intensity. The class is
    initialized from a class instance of SHMagCoeffs using the method
    expand().

    Attributes:

    rad            : SHGrid class instance of the radial component of the
                     magnetic field evaluated on an ellipsoid.
    theta          : SHGrid class instance of the theta component of the
                     magnetic field evaluated on an ellipsoid.
    phi            : SHGrid class instance of the phi component of the
                     magnetic field evaluated on an ellipsoid.
    total          : SHGrid class instance of the total magnetic intensity on
                     an ellipsoid.
    pot            : SHGrid class instance of the magnetic potential
                     evaluated on an ellipsoid.
    a              : Semimajor axis of the reference ellipsoid.
    f              : Flattening of the reference ellipsoid, f=(a-b)/a.
    lmax           : The maximum spherical harmonic degree resolvable by the
                     grids.
    lmax_calc      : The maximum spherical harmonic degree of the magnetic
                     potential used in creating the grids.
    nlat, nlon     : The number of latitude and longitude bands in the grids.
    sampling       : The longitudinal sampling scheme of the grids: either
                     1 for nlon=nlat or 2 for nlon=2*nlat.

    Methods:

    plot()        : Plot all three components of the magnetic field and the
                    total magnetic intensity.
    plot_rad()    : Plot the radial component of the magnetic field.
    plot_theta()  : Plot the theta component of the magnetic field.
    plot_phi()    : Plot the phi component of the magnetic field.
    plot_total()  : Plot the total magnetic intensity.
    plot_pot()    : Plot the magnetic potential.
    copy()        : Return a copy of the class instance.
    info()        : Print a summary of the data stored in the SHMagGrid
                    instance.
    """

    def __init__(self, rad, theta, phi, total, pot, a, f, lmax, lmax_calc):
        """
        Initialize the SHMagGrid class.
        """
        self.rad = _SHGrid.from_array(rad, grid='DH')
        self.theta = _SHGrid.from_array(theta, grid='DH')
        self.phi = _SHGrid.from_array(phi, grid='DH')
        self.total = _SHGrid.from_array(total, grid='DH')
        self.pot = _SHGrid.from_array(pot, grid='DH')
        self.grid = self.rad.grid
        self.sampling = self.rad.sampling
        self.nlat = self.rad.nlat
        self.nlon = self.rad.nlon
        self.a = a
        self.f = f
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
        Print a summary of the data stored in the SHMagGrid class instance.

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
                'lmax = {:d}\n'
                'lmax_calc = {:d}\n'
                'a (m)= {:e}\n'
                'f = {:e}'
                .format(self.nlat, self.nlon, self.lmax, self.lmax_calc,
                        self.a, self.f))
        return str

    def plot_rad(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$B_r$, nT', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the radial component of the magnetic field.

        Usage
        -----
        x.plot_rad([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname, **kwargs])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = '$B_r$, nT'
            Text label for the colorbar.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        kwargs : optional
            Keyword arguements that will be sent to the SHGrid.plot()
            and plt.imshow() methods.
        """
        if ax is None:
            fig, axes = self.rad.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, show=False, **kwargs)
            if show:
                fig.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.rad.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_theta(self, colorbar=True, cb_orientation='vertical',
                   cb_label='$B_\\theta$, nT', ax=None, show=True,
                   fname=None, **kwargs):
        """
        Plot the theta component of the magnetic field.

        Usage
        -----
        x.plot_theta([tick_interval, xlabel, ylabel, ax, colorbar,
                      cb_orientation, cb_label, show, fname, **kwargs])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = '$B_\\theta$, nT'
            Text label for the colorbar.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        kwargs : optional
            Keyword arguements that will be sent to the SHGrid.plot()
            and plt.imshow() methods.
        """
        if ax is None:
            fig, axes = self.theta.plot(colorbar=colorbar,
                                        cb_orientation=cb_orientation,
                                        cb_label=cb_label, show=False,
                                        **kwargs)
            if show:
                fig.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.theta.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                            cb_label=cb_label, ax=ax, **kwargs)

    def plot_phi(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$B_\phi$, nT', ax=None, show=True,
                 fname=None, **kwargs):
        """
        Plot the phi component of the magnetic field.

        Usage
        -----
        x.plot_phi([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname, **kwargs])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = '$B_\phi$, nT'
            Text label for the colorbar.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        kwargs : optional
            Keyword arguements that will be sent to the SHGrid.plot()
            and plt.imshow() methods.
        """
        if ax is None:
            fig, axes = self.phi.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, show=False, **kwargs)
            if show:
                fig.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.phi.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_total(self, colorbar=True, cb_orientation='vertical',
                   cb_label='$|B|$, nT', ax=None, show=True, fname=None,
                   **kwargs):
        """
        Plot the total magnetic intensity.

        Usage
        -----
        x.plot_total([tick_interval, xlabel, ylabel, ax, colorbar,
                      cb_orientation, cb_label, show, fname, **kwargs])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = '$|B|$, nT'
            Text label for the colorbar.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        kwargs : optional
            Keyword arguements that will be sent to the SHGrid.plot()
            and plt.imshow() methods.
        """
        if ax is None:
            fig, axes = self.total.plot(
                colorbar=colorbar, cb_orientation=cb_orientation,
                cb_label=cb_label, show=False, **kwargs)

            if show:
                fig.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.total.plot(
                colorbar=colorbar, cb_orientation=cb_orientation,
                cb_label=cb_label, ax=ax, **kwargs)

    def plot_pot(self, colorbar=True, cb_orientation='vertical',
                 cb_label='Potential, nT m', ax=None, show=True,
                 fname=None, **kwargs):
        """
        Plot the magnetic potential.

        Usage
        -----
        x.plot_pot([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname, **kwargs])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = [30, 30]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = 'potential, nT m'
            Text label for the colorbar.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        kwargs : optional
            Keyword arguements that will be sent to the SHGrid.plot()
            and plt.imshow() methods.
        """
        if ax is None:
            fig, axes = self.pot.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, show=False, **kwargs)
            if show:
                fig.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.pot.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot(self, colorbar=True, cb_orientation='horizontal',
             tick_interval=[60, 60], minor_tick_interval=[20, 20],
             xlabel='Longitude', ylabel='Latitude',
             axes_labelsize=9, tick_labelsize=8, show=True, fname=None,
             **kwargs):
        """
        Plot the three vector components of the magnetic field and the total
        magnetic intensity.

        Usage
        -----
        x.plot([tick_interval, minor_tick_interval, xlabel, ylabel,
                colorbar, cb_orientation, cb_label, axes_labelsize,
                tick_labelsize, show, fname, **kwargs])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = [60, 60]
            Intervals to use when plotting the major x and y ticks. If set to
            None, major ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [20, 20]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        xlabel : str, optional, default = 'Longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'Latitude'
            Label for the latitude axis.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = None
            Text label for the colorbar.
        axes_labelsize : int, optional, default = 9
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = 8
            The font size for the x and y tick labels.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        kwargs : optional
            Keyword arguements that will be sent to the SHGrid.plot()
            and plt.imshow() methods.
        """
        if colorbar is True:
            if cb_orientation == 'horizontal':
                scale = 0.8
            else:
                scale = 0.5
        else:
            scale = 0.6
        figsize = (_mpl.rcParams['figure.figsize'][0],
                    _mpl.rcParams['figure.figsize'][0] * scale)

        fig, ax = _plt.subplots(2, 2, figsize=figsize)
        self.plot_rad(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[0], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel,
                      axes_labelsize=axes_labelsize,
                      tick_labelsize=tick_labelsize,
                      minor_tick_interval=minor_tick_interval,
                      **kwargs)
        self.plot_theta(colorbar=colorbar, cb_orientation=cb_orientation,
                        ax=ax.flat[1], tick_interval=tick_interval,
                        xlabel=xlabel, ylabel=ylabel,
                        axes_labelsize=axes_labelsize,
                        tick_labelsize=tick_labelsize,
                        minor_tick_interval=minor_tick_interval,
                        **kwargs)
        self.plot_phi(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[2], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel,
                      axes_labelsize=axes_labelsize,
                      minor_tick_interval=minor_tick_interval,
                      tick_labelsize=tick_labelsize,**kwargs)
        self.plot_total(colorbar=colorbar, cb_orientation=cb_orientation,
                        ax=ax.flat[3], tick_interval=tick_interval,
                        xlabel=xlabel, ylabel=ylabel,
                        axes_labelsize=axes_labelsize,
                        tick_labelsize=tick_labelsize,
                        minor_tick_interval=minor_tick_interval,
                        **kwargs)
        fig.tight_layout(pad=0.5)

        if show:
            fig.show()

        if fname is not None:
            fig.savefig(fname)
        return fig, ax
