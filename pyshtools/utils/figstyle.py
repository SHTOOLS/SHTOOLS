"""
Set matplotlib parameters for creating publication quality graphics.
"""
import matplotlib.pyplot as _plt


def figstyle(rel_width=0.75, screen_dpi=114, aspect_ratio=4/3,
             max_width=7.48031):
    """
    Set matplotlib parameters for creating publication quality graphics.

    Usage
    -----
    figstyle([rel_width, screen_dpi, aspect_ratio, max_width])

    Parameters
    ----------
    rel_width : float, optional, default = 0.75
        The relative width of the plot (from 0 to 1) wih respect to max_width.
    screen_dpi : int, optional, default = 114
        The screen resolution of the display in dpi, which determines the
        size of the plot on the display.
    aspect_ratio : float, optional, default = 4/3
        The aspect ratio of the plot.
    max_width : float, optional, default = 7.48031
        The maximum width of the usable area of a journal page in inches.

    Description
    -----------
    This function sets a variety of matplotlib parameters for creating
    publication quality graphics. The default parameters are tailored to
    AGU/Wiley-Blackwell journals that accept relative widths of 0.5, 0.75,
    and 1. To reset the maplotlib parameters to their default values, use

        matplotlib.pyplot.style.use('default')
    """
    width_x = max_width * rel_width
    width_y = max_width * rel_width / aspect_ratio

    shtools = {
        # fonts
        'font.size': 10,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Myriad Pro', 'DejaVu Sans',
                            'Bitstream Vera Sans',
                            'Verdana', 'Arial', 'Helvetica'],
        'axes.titlesize': 10,
        'axes.labelsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 9,
        'text.usetex': False,
        'axes.formatter.limits': (-3, 3),
        # figure
        'figure.dpi': screen_dpi,
        'figure.figsize': (width_x, width_y),
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
        # images
        'image.lut': 65536,  # 16 bit
        # savefig
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'savefig.dpi': 600,
        'savefig.format': 'pdf'
        }

    _plt.style.use([shtools])
