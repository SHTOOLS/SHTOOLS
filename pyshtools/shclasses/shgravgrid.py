"""
    Class for grids of the three components of the graivty field, the
    gravitational disturbance, and the gravitational potential.
"""
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy

from .shcoeffsgrid import SHGrid as _SHGrid


class SHGravGrid(object):
    """
    Class for grids of the gravitational potential, three components of the
    gravity field, and the total gravitational disturbance. The class is
    initialized using SHGravCoeffs.expand().

    Attributes

    rad            : SHGrid class instance of the radial component of the
                     gravitational acceleration evaluated on an ellipsoid.
    theta          : SHGrid class instance of the theta component of the
                     gravitational acceleration evaluated on an ellipsoid..
    phi            : SHGrid class instance of the phi component of the
                     gravitational acceleration evaluated on an ellipsoid..
    total          : SHGrid class instance of the total gravitational
                     acceleration with the normal gravity removed, on an
                     ellipsoid.
    pot            : SHGrid class instance of the gravitational potential
                     evaluated on an ellipsoid.
    gm             : Gravitational constant time the mass of the body.
    a              : Semimajor axis of the reference ellipsoid.
    f              : Flattening of the reference ellipsoid,
                     f=(R_equator-R_pole)/R_equator.
    omega          : Angular rotation rate of the body.
    normal_gravity : True if the normal gravity is removed from the total
                     gravitational acceleration.
    lmax           : The maximum spherical harmonic degree resolvable by the
                     grids.
    lmax_calc      : The maximum spherical harmonic degree of the gravitational
                     potential used in creating the grids.
    nlat, nlon     : The number of latitude and longitude bands in the grids.
    sampling       : The longitudinal sampling scheme of the grids: either 1
                     for nlong=nlat or 2 for nlong=2*nlat.

    Methods

    plot()        : Plot all three components of the gravity field with the
                    total gravity disturbance.
    plot_rad()    : Plot the radial component of the gravity field.
    plot_theta()  : Plot the theta component of the gravity field.
    plot_phi()    : Plot the phi component of the gravity field.
    plot_total()  : Plot the total gravity disturbance.
    plot_pot()    : Plot the gravitational potential.
    copy()        : Return a copy of the class instance.
    info()        : Print a summary of the data stored in the SHGrid instance.
    """

    def __init__(self, rad, theta, phi, total, pot, gm, a, f, omega,
                 normal_gravity, lmax, lmax_calc):
        """
        Initialize the SHGravGradGrid class.
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
        self.gm = gm
        self.a = a
        self.f = f
        self.omega = omega
        self.normal_gravity = normal_gravity
        self.lmax = lmax
        self.lmax_calc = lmax_calc

    def copy(self):
        """Return a deep copy of the class instance."""
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
        str = ('grid = {:s}\n'.format(repr(self.grid)))
        if self.grid == 'DH':
            str += 'sampling = {:d}\n'.format(self.sampling)
        str += ('nlat = {:d}\n'
                'nlon = {:d}\n'
                'lmax = {:d}\n'
                'lmax_calc = {:d}\n'
                'gm (m3 / s2) = {:e}\n'
                'a (m)= {:e}\n'
                'f = {:e}\n'
                'omega (rad / s) = {:s}\n'
                'normal gravity removed = {:s}'
                .format(self.nlat, self.nlon, self.lmax, self.lmax_calc,
                        self.gm, self.a, self.f, repr(self.omega),
                        repr(self.normal_gravity)))
        return str

    def plot_rad(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$g_r$,  m s$^{-2}$', **kwargs):
        """
        Plot the radial component of the gravity field.

        Usage
        -----
        x.plot_rad([tick_interval, ax, colorbar, cb_orientation, cb_label,
                    show, fname])

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
        colorbar : bool, optional, default = False
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar; either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = None
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
        self.rad.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                      cb_label=cb_label, **kwargs)

    def plot_theta(self, colorbar=True, cb_orientation='vertical',
                   cb_label='$g_\\theta$,  m s$^{-2}$', **kwargs):
        """
        Plot the theta component of the gravity field.

        Usage
        -----
        x.plot_theta([tick_interval, ax, colorbar, cb_orientation, cb_label,
                      show, fname])

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
        colorbar : bool, optional, default = False
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar; either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = None
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
        self.theta.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                        cb_label=cb_label, **kwargs)

    def plot_phi(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$g_\phi$,  m s$^{-2}$', **kwargs):
        """
        Plot the phi component of the gravity field.

        Usage
        -----
        x.plot_phi([tick_interval, ax, colorbar, cb_orientation, cb_label,
                    show, fname])

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
        colorbar : bool, optional, default = False
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar; either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = None
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
        self.phi.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                      cb_label=cb_label, **kwargs)

    def plot_total(self, colorbar=True, cb_orientation='vertical',
                   cb_label='gravity disturbance, m s$^{-2}$', **kwargs):
        """
        Plot the total gravity disturbance.

        Usage
        -----
        x.plot_total([tick_interval, ax, colorbar, cb_orientation, cb_label,
                      show, fname])

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
        colorbar : bool, optional, default = False
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'horizontal'
            Orientation of the colorbar; either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = None
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
        if self.normal_gravity is True:
            cb_label = 'gravity disturbance, mGal'
            (self.total*1.e5).plot(colorbar=colorbar,
                                   cb_orientation=cb_orientation,
                                   cb_label=cb_label, **kwargs)
        else:
            self.total.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                            cb_label=cb_label, **kwargs)

    def plot_pot(self, colorbar=True, cb_orientation='vertical',
                 cb_label='potential,  m s$^{-1}$', **kwargs):
        """
        Plot the gravitational potential.

        Usage
        -----
        x.plot_pot([tick_interval, ax, colorbar, cb_orientation, cb_label,
                    show, fname])

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
        colorbar : bool, optional, default = False
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar; either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = None
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
        self.rad.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                      cb_label=cb_label, **kwargs)

    def plot(self, colorbar=True, cb_orientation='horizontal',
             tick_interval=[45, 45], xlabel='', ylabel='', show=True,
             fname=None):
        """
        Plot the three components of the gravity field and the gravity
        disturbance.

        Usage
        -----
        x.plot([tick_interval, ax, colorbar, cb_orientation, cb_label,
                show, fname])

        tick_interval : list or tuple, optional, default = [45, 45]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
            Label for the latitude axis.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        colorbar : bool, optional, default = False
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar; either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = None
            Text label for the colorbar.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        """
        fig, ax = _plt.subplots(2, 2)
        self.plot_rad(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[0], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel)
        self.plot_theta(colorbar=colorbar, cb_orientation=cb_orientation,
                        ax=ax.flat[1], tick_interval=tick_interval,
                        xlabel=xlabel, ylabel=ylabel)
        self.plot_phi(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[2], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel)
        self.plot_total(colorbar=colorbar, cb_orientation=cb_orientation,
                        ax=ax.flat[3], tick_interval=tick_interval,
                        xlabel=xlabel, ylabel=ylabel)
        fig.tight_layout(pad=0.5)

        if ax is None:
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, ax