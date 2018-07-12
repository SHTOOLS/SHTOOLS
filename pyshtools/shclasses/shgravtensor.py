"""
    Class for the gravity 'gradient' tensor.
"""
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy
from scipy.linalg import eigvalsh as _eigvalsh

from .shcoeffsgrid import SHGrid as _SHGrid


class SHGravTensor(object):
    """
    Class for the gravity 'gradient' tensor and eigenvalues. The class is
    initialized from a class instance of SHGravCoeffs using the method
    tensor().

    Attributes:

    vxx, vxy, vzz,   : The 9 components of the gravity tensor.
    vyx, vyy, vyz,
    vzx, vzy, vzz
    i0, i1, i2, i    : The three invariants of the gravity tensor and a
                       derived quantity that is bounded between 0 and 1.
    gm               : The gravitational constant times the mass of the body.
    a                : Semimajor axis of the reference ellipsoid.
    f                : Flattening of the reference ellipsoid, f=(a-b)/a.
    lmax             : The maximum spherical harmonic degree resolvable by the
                       grids.
    lmax_calc        : The maximum spherical harmonic degree of the
                       gravitational potential used in creating the grids.
    nlat, nlon       : The number of latitude and longitude bands in the grids.
    sampling         : The longitudinal sampling scheme of the grids: either 1
                       for nlong=nlat or 2 for nlong=2*nlat.

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

    copy()         : Return a copy of the class instance.
    info()         : Print a summary of the data stored in the SHGravTensor
                     instance.
    """

    def __init__(self, vxx, vyy, vzz, vxy, vxz, vyz, gm, a, f, lmax,
                 lmax_calc):
        """
        Initialize the SHGravTensor class.
        """
        self.vxx = _SHGrid.from_array(vxx, grid='DH')
        self.vyy = _SHGrid.from_array(vyy, grid='DH')
        self.vzz = _SHGrid.from_array(vzz, grid='DH')
        self.vxy = _SHGrid.from_array(vxy, grid='DH')
        self.vxz = _SHGrid.from_array(vxz, grid='DH')
        self.vyz = _SHGrid.from_array(vyz, grid='DH')
        self.vyx = self.vxy
        self.vzx = self.vxz
        self.vzy = self.vyz
        self.grid = self.vxx.grid
        self.sampling = self.vxx.sampling
        self.nlat = self.vxx.nlat
        self.nlon = self.vxx.nlon
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

    def compute_invar(self):
        """
        Compute the three invariants (I0, I1, I2) of the gravity tensor, as
        well as the quantity I = -(I2/2)**2 / (I1/3)**3.
        """
        self.i0 = self.vxx + self.vyy + self.vzz
        self.i1 = (self.vxx*self.vyy + self.vyy*self.vzz + self.vxx*self.vzz -
                   self.vxy**2 - self.vyz**2 - self.vxz**2)
        self.i2 = (self.vxx*(self.vyy*self.vzz - self.vyz**2) +
                   self.vxy*(self.vyz*self.vxz - self.vxy*self.vzz) +
                   self.vxz*(self.vxy*self.vyz - self.vxz*self.vyy))
        self.i = (-1.) * (self.i2 / 2.)**2
        self.i.data[1:, :] /= (self.i1.data[1:, :] / 3.)**3

    def compute_eig(self):
        """
        Compute the three eigenvalues of the gravity tensor: eig1, eig2, ei3.
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
        Compute the two horizontal eigenvalues of the gravity tensor (eigh1,
        and eigh2), as well as the combined maximum absolute value of the
        two (eighh).
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
                'f = {:e}'
                .format(self.nlat, self.nlon, self.lmax, self.lmax_calc,
                        self.gm, self.a, self.f))
        return str

    def plot_vxx(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{xx}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vxx component of the gravity tensor.

        Usage
        -----
        x.plot_vxx([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
        cb_label : str, optional, default = '$V_{xx}$, Eotvos'
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
            fig, axes = self.vxx.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vxx.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vyy(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{yy}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vyy component of the gravity tensor.

        Usage
        -----
        x.plot_vyy([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{yy}$, Eotvos'
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
            fig, axes = self.vyy.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vyy.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vzz(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{zz}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vzz component of the gravity tensor.

        Usage
        -----
        x.plot_vzz([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{zz}$, Eotvos'
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
            fig, axes = self.vzz.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vzz.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vxy(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{xy}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vxy component of the gravity tensor.

        Usage
        -----
        x.plot_vxy([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{xy}$, Eotvos'
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
            fig, axes = self.vxy.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vxy.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vyx(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{yx}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vyx component of the gravity tensor.

        Usage
        -----
        x.plot_vyx([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{yx}$, Eotvos'
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
            fig, axes = self.vyx.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vyx.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vxz(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{xz}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vxz component of the gravity tensor.

        Usage
        -----
        x.plot_vxz([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{xz}$, Eotvos'
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
            fig, axes = self.vxz.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vxz.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vzx(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{zx}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vzx component of the gravity tensor.

        Usage
        -----
        x.plot_vzx([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{zx}$, Eotvos'
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
            fig, axes = self.vzx.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vzx.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vyz(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{yz}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vxx component of the gravity tensor.

        Usage
        -----
        x.plot_vyz([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{yz}$, Eotvos'
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
            fig, axes = self.vyz.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vyz.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot_vzy(self, colorbar=True, cb_orientation='vertical',
                 cb_label='$V_{xx}$, Eotvos', ax=None, show=True, fname=None,
                 **kwargs):
        """
        Plot the Vzy component of the gravity tensor.

        Usage
        -----
        x.plot_vzy([tick_interval, xlabel, ylabel, ax, colorbar,
                    cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$V_{zy}$, Eotvos'
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
            fig, axes = self.vzy.plot(colorbar=colorbar,
                                      cb_orientation=cb_orientation,
                                      cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.vzy.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                          cb_label=cb_label, ax=ax, **kwargs)

    def plot(self, colorbar=True, cb_orientation='horizontal',
             tick_interval=None, xlabel='', ylabel='', show=True,
             fname=None, **kwargs):
        """
        Plot the 9 components of the gravity 'gradient' tensor.

        Usage
        -----
        x.plot([tick_interval, xlabel, ylabel, ax, colorbar, cb_orientation,
                    cb_label, show, fname])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = None
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
            Label for the latitude axis.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'horizontal'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
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
        fig, ax = _plt.subplots(3, 3)
        self.plot_vxx(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[0], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vxy(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[1], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vxz(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[2], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vyx(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[3], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vyy(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[4], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vyz(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[5], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vzx(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[6], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vzy(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[7], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_vzz(colorbar=colorbar, cb_orientation=cb_orientation,
                      ax=ax.flat[8], tick_interval=tick_interval,
                      xlabel=xlabel, ylabel=ylabel, **kwargs)

        fig.tight_layout(pad=0.5)

        if show:
            _plt.show()

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def plot_i0(self, colorbar=True, cb_orientation='vertical',
                cb_label='Tr $V_{ij}$, Eotvos', ax=None, show=True,
                fname=None, **kwargs):
        """
        Plot the first invariant I0 (the trace) of the gravity tensor

            I0 = vxx + vyy + vzz

        which should be identically zero.

        Usage
        -----
        x.plot_i0([tick_interval, xlabel, ylabel, ax, colorbar, cb_orientation,
                   cb_label, show, fname])

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
        cb_label : str, optional, default = 'Tr $V_{ij}$, Eotvos'
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
        if self.i0 is None:
            self.compute_invar()

        if ax is None:
            fig, axes = self.i0.plot(colorbar=colorbar,
                                     cb_orientation=cb_orientation,
                                     cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.i0.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                         cb_label=cb_label, ax=ax, **kwargs)

    def plot_i1(self, colorbar=True, cb_orientation='vertical',
                cb_label='$I_1$, Eotvos$^2$', ax=None, show=True, fname=None,
                **kwargs):
        """
        Plot the second invariant I1 of the gravity tensor:

            I1 = vxx*vyy + vyy*vzz + vxx*vzz - vxy**2 - vyz**2 - vxz**2

        Usage
        -----
        x.plot_i1([tick_interval, xlabel, ylabel, ax, colorbar, cb_orientation,
                   cb_label, show, fname])

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
        cb_label : str, optional, default = '$I_1$, Eotvos$^2$'
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
        if self.i1 is None:
            self.compute_invar()

        if ax is None:
            fig, axes = self.i1.plot(colorbar=colorbar,
                                     cb_orientation=cb_orientation,
                                     cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.i1.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                         cb_label=cb_label, ax=ax, **kwargs)

    def plot_i2(self, colorbar=True, cb_orientation='vertical',
                cb_label='det $V_{ij}$, Eotvos$^3$', ax=None, show=True,
                fname=None, **kwargs):
        """
        Plot the third invariant I2 (the determinant) of the gravity tensor:

           I2 = vxx*(vyy*vzz - vyz**2) + vxy*(vyz*vxz - vxy*vzz)
                + vxz*(vxy*vyz - vxz*vyy)

        Usage
        -----
        x.plot_i2([tick_interval, xlabel, ylabel, ax, colorbar, cb_orientation,
                   cb_label, show, fname])

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
        cb_label : str, optional, default = 'det $V_{ij}$, Eotvos$^3$'
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
        if self.i2 is None:
            self.compute_invar()

        if ax is None:
            fig, axes = self.i2.plot(colorbar=colorbar,
                                     cb_orientation=cb_orientation,
                                     cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.i2.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                         cb_label=cb_label, ax=ax, **kwargs)

    def plot_i(self, colorbar=True, cb_orientation='vertical',
               cb_label='$-(I_2/2)^{2} / (I_1/3)^{3}$, Eotvos$^{-1}$',
               ax=None, show=True, fname=None, **kwargs):
        """
        Plot the quantity I of Pedersen and Rasmussen (1990)

           I = -(I2/2)**2 / (I1/3)**3

        that is bounded by 0 and 1.

        Usage
        -----
        x.plot_i([tick_interval, xlabel, ylabel, ax, colorbar, cb_orientation,
                  cb_label, show, fname])

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
        cb_label : str, optional, default = '$-(I_2/2)^{2} / (I_1/3)^{3}$,
                                             Eotvos$^{-1}$'
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
        if self.i is None:
            self.compute_invar()

        if ax is None:
            fig, axes = self.i.plot(colorbar=colorbar,
                                    cb_orientation=cb_orientation,
                                    cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.i.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                        cb_label=cb_label, ax=ax, **kwargs)

    def plot_invar(self, colorbar=True, cb_orientation='horizontal',
                   tick_interval=None, xlabel='', ylabel='', show=True,
                   fname=None, **kwargs):
        """
        Plot the three invariants of the gravity tensor and the derived
        quantity I.

        Usage
        -----
        x.plot_invar([tick_interval, xlabel, ylabel, ax, colorbar,
                      cb_orientation, cb_label, show, fname])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = None
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
            Label for the latitude axis.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'horizontal'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
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
        fig, ax = _plt.subplots(2, 2)
        self.plot_i0(colorbar=colorbar, cb_orientation=cb_orientation,
                     ax=ax.flat[0], tick_interval=tick_interval,
                     xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_i1(colorbar=colorbar, cb_orientation=cb_orientation,
                     ax=ax.flat[1], tick_interval=tick_interval,
                     xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_i2(colorbar=colorbar, cb_orientation=cb_orientation,
                     ax=ax.flat[2], tick_interval=tick_interval,
                     xlabel=xlabel, ylabel=ylabel, **kwargs)
        self.plot_i(colorbar=colorbar, cb_orientation=cb_orientation,
                    ax=ax.flat[3], tick_interval=tick_interval,
                    xlabel=xlabel, ylabel=ylabel, **kwargs)

        fig.tight_layout(pad=0.5)

        if show:
            _plt.show()

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def plot_eig1(self, colorbar=True, cb_orientation='vertical',
                  cb_label='$\lambda_1$, Eotvos$^{-1}$', ax=None, show=True,
                  fname=None, **kwargs):
        """
        Plot the first eigenvalue of the gravity tensor.

        Usage
        -----
        x.plot_eig1([tick_interval, xlabel, ylabel, ax, colorbar,
                     cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$\lambda_1$, Eotvos$^{-1}$'
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
        if self.eig1 is None:
            self.compute_eig()

        if ax is None:
            fig, axes = self.eig1.plot(colorbar=colorbar,
                                       cb_orientation=cb_orientation,
                                       cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.eig1.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                           cb_label=cb_label, ax=ax, **kwargs)

    def plot_eig2(self, colorbar=True, cb_orientation='vertical',
                  cb_label='$\lambda_2$, Eotvos$^{-1}$', ax=None, show=True,
                  fname=None, **kwargs):
        """
        Plot the second eigenvalue of the gravity tensor.

        Usage
        -----
        x.plot_eig2([tick_interval, xlabel, ylabel, ax, colorbar,
                     cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$\lambda_2$, Eotvos$^{-1}$'
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
        if self.eig2 is None:
            self.compute_eig()

        if ax is None:
            fig, axes = self.eig2.plot(colorbar=colorbar,
                                       cb_orientation=cb_orientation,
                                       cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.eig2.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                           cb_label=cb_label, ax=ax, **kwargs)

    def plot_eig3(self, colorbar=True, cb_orientation='vertical',
                  cb_label='$\lambda_3$, Eotvos$^{-1}$', ax=None, show=True,
                  fname=None, **kwargs):
        """
        Plot the third eigenvalue of the gravity tensor.

        Usage
        -----
        x.plot_eig3([tick_interval, xlabel, ylabel, ax, colorbar,
                     cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$\lambda_3$, Eotvos$^{-1}$'
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
        if self.eig3 is None:
            self.compute_eig()

        if ax is None:
            fig, axes = self.eig3.plot(colorbar=colorbar,
                                       cb_orientation=cb_orientation,
                                       cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.eig3.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                           cb_label=cb_label, ax=ax, **kwargs)

    def plot_eigs(self, colorbar=True, cb_orientation='vertical',
                  tick_interval=[45, 45], xlabel='', ylabel='', show=True,
                  fname=None, **kwargs):
        """
        Plot the three eigenvalues of the gravity tensor.

        Usage
        -----
        x.plot_eigs([tick_interval, xlabel, ylabel, ax, colorbar,
                     cb_orientation, cb_label, show, fname])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = None
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
            Label for the latitude axis.
        colorbar : bool, optional, default = False
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
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
        fig, ax = _plt.subplots(3, 1)

        self.plot_eig1(colorbar=colorbar, cb_orientation=cb_orientation,
                       ax=ax.flat[0], xlabel=xlabel, ylabel=ylabel,
                       tick_interval=tick_interval, **kwargs)
        self.plot_eig2(colorbar=colorbar, cb_orientation=cb_orientation,
                       ax=ax.flat[1], xlabel=xlabel, ylabel=ylabel,
                       tick_interval=tick_interval, **kwargs)
        self.plot_eig3(colorbar=colorbar, cb_orientation=cb_orientation,
                       ax=ax.flat[2], xlabel=xlabel, ylabel=ylabel,
                       tick_interval=tick_interval, **kwargs)

        fig.tight_layout(pad=0.5)

        if show:
            _plt.show()

        if fname is not None:
            fig.savefig(fname)
        return fig, ax

    def plot_eigh1(self, colorbar=True, cb_orientation='vertical',
                   cb_label='$\lambda_{h1}$, Eotvos$^{-1}$', ax=None,
                   show=True, fname=None, **kwargs):
        """
        Plot the first eigenvalue of the horizontal gravity tensor.

        Usage
        -----
        x.plot_eigh1([tick_interval, xlabel, ylabel, ax, colorbar,
                      cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$\lambda_{h1}$, Eotvos$^{-1}$'
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
        if self.eigh1 is None:
            self.compute_eigh()

        if ax is None:
            fig, axes = self.eigh1.plot(colorbar=colorbar,
                                        cb_orientation=cb_orientation,
                                        cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.eigh1.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                            cb_label=cb_label, ax=ax, **kwargs)

    def plot_eigh2(self, colorbar=True, cb_orientation='vertical',
                   cb_label='$\lambda_{h2}$, Eotvos$^{-1}$', ax=None,
                   show=True, fname=None, **kwargs):
        """
        Plot the second eigenvalue of the horizontal gravity tensor.

        Usage
        -----
        x.plot_eigh2([tick_interval, xlabel, ylabel, ax, colorbar,
                      cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$\lambda_{h2}$, Eotvos$^{-1}$'
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
        if self.eigh2 is None:
            self.compute_eigh()

        if ax is None:
            fig, axes = self.eigh2.plot(colorbar=colorbar,
                                        cb_orientation=cb_orientation,
                                        cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.eigh2.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                            cb_label=cb_label, ax=ax, **kwargs)

    def plot_eighh(self, colorbar=True, cb_orientation='vertical',
                   cb_label='$\lambda_{hh}$, Eotvos$^{-1}$', ax=None,
                   show=True, fname=None, **kwargs):
        """
        Plot the maximum absolute value eigenvalue of the horizontal gravity
        tensor.

        Usage
        -----
        x.plot_eighh([tick_interval, xlabel, ylabel, ax, colorbar,
                      cb_orientation, cb_label, show, fname])

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
        cb_label : str, optional, default = '$\lambda_{hh}$, Eotvos$^{-1}$'
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
        if self.eighh is None:
            self.compute_eigh()

        if ax is None:
            fig, axes = self.eighh.plot(colorbar=colorbar,
                                        cb_orientation=cb_orientation,
                                        cb_label=cb_label, **kwargs)
            if show:
                _plt.show()

            if fname is not None:
                fig.savefig(fname)
            return fig, axes

        else:
            self.eighh.plot(colorbar=colorbar, cb_orientation=cb_orientation,
                            cb_label=cb_label, ax=ax, **kwargs)

    def plot_eigh(self, colorbar=True, cb_orientation='vertical',
                  tick_interval=[45, 45], xlabel='', ylabel='', show=True,
                  fname=None, **kwargs):
        """
        Plot the two eigenvalues and maximum absolute value eigevnalue of the
        horizontal gravity tensor.

        Usage
        -----
        x.plot_eigs([tick_interval, xlabel, ylabel, ax, colorbar,
                     cb_orientation, cb_label, show, fname])

        Parameters
        ----------
        tick_interval : list or tuple, optional, default = None
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = ''
            Label for the longitude axis.
        ylabel : str, optional, default = ''
            Label for the latitude axis.
        colorbar : bool, optional, default = True
            If True, plot a colorbar.
        cb_orientation : str, optional, default = 'vertical'
            Orientation of the colorbar: either 'vertical' or 'horizontal'.
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
        fig, ax = _plt.subplots(3, 1)

        self.plot_eigh1(colorbar=colorbar, cb_orientation=cb_orientation,
                        ax=ax.flat[0], xlabel=xlabel, ylabel=ylabel,
                        tick_interval=tick_interval, **kwargs)
        self.plot_eigh2(colorbar=colorbar, cb_orientation=cb_orientation,
                        ax=ax.flat[1], xlabel=xlabel, ylabel=ylabel,
                        tick_interval=tick_interval, **kwargs)
        self.plot_eighh(colorbar=colorbar, cb_orientation=cb_orientation,
                        ax=ax.flat[2], xlabel=xlabel, ylabel=ylabel,
                        tick_interval=tick_interval, **kwargs)

        fig.tight_layout(pad=0.5)

        if show:
            _plt.show()

        if fname is not None:
            fig.savefig(fname)
        return fig, ax
