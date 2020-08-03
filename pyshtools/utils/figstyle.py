"""
Set matplotlib parameters for creating publication quality graphics.

Note: SHGrid.plotgmt() sets the title, axes, and tick font sizes using the
values defined by mpl.rcParams (which are modified when utils.figstyle is
called).
"""
import matplotlib.pyplot as _plt


def figstyle(rel_width=0.75, screen_dpi=114, aspect_ratio=4/3,
             max_width=7.48031, figsize=None, units='i'):
    """
    Set matplotlib parameters for creating publication quality graphics.

    Usage
    -----
    figstyle([rel_width, screen_dpi, aspect_ratio, max_width, figsize, units])

    Parameters
    ----------
    rel_width : float, optional, default = 0.75
        The relative width of the plot (from 0 to 1) wih respect to width.
    screen_dpi : int, optional, default = 114
        The screen resolution of the display in dpi, which determines the
        size of the plot on the display.
    aspect_ratio : float, optional, default = 4/3
        The aspect ratio (horizontal/vertical) of the plot.
    max_width : float, optional, default = 7.48031
        The maximum width of the usable area of a journal page in inches.
    figsize : list, optional, default = None
        A list containing the width and height of the figure.
    units : str, optional, default = 'i'
        The measurement units of max_width and figsize, either 'i' for inches
        or 'c' or 'cm' for centimeters.

    Notes
    -----
    This function sets a variety of matplotlib parameters for creating
    publication quality graphics. The default parameters are tailored to
    AGU/Wiley-Blackwell journals that accept relative widths of 0.5, 0.75,
    and 1. To reset the maplotlib parameters to their default values, use

        matplotlib.pyplot.style.use('default')
    """
    if figsize is None:
        figsize = (max_width * rel_width, max_width * rel_width / aspect_ratio)

    if units == 'c' or units == 'cm':
        # convert to inches
        figsize /= 2.54

    shtools = {
        # fonts
        'font.size': 10,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Myriad Pro', 'DejaVu Sans',
                            'Bitstream Vera Sans',
                            'Verdana', 'Arial', 'Helvetica'],
        'axes.titlesize': 11,
        'axes.titlepad': 10,
        'axes.labelsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 9,
        'text.usetex': False,
        'axes.formatter.limits': (-3, 4),
        'mathtext.fontset': 'stix',
        # figure
        'figure.dpi': screen_dpi,
        'figure.figsize': figsize,
        # line and tick widths
        'axes.linewidth': 1,
        'lines.linewidth': 1.5,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.minor.width': 0.6,
        'xtick.minor.width': 0.6,
        'xtick.top': True,
        'ytick.right': True,
        # grids
        'grid.linewidth': 0.3,
        'grid.color': 'k',
        'grid.linestyle': '-',
        # legends
        'legend.framealpha': 1.,
        'legend.edgecolor': 'k',
        # error bars
        'errorbar.capsize': 1.5,
        # images
        'image.lut': 65536,  # 16 bit
        # savefig
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'savefig.dpi': 600,
        'savefig.format': 'pdf'
        }

    _plt.style.use([shtools])
