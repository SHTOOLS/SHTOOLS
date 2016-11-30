"""
Contains configuration parameters for matplotlib.
"""

from __future__ import absolute_import, division, print_function

# ---- some typical sizes for an A4 text document ----
pt2inch = 1.0 / 72.27
textwidth_pt = 504.0
columnwidth_pt = 240.0
textwidth_in = textwidth_pt * pt2inch
columnwidth_in = columnwidth_pt * pt2inch

# ---- this is the main configuration dictionary ----
style_shtools = {'font.size': 7,
                 'axes.labelsize': 'medium',
                 'axes.titlesize': 'medium',
                 'xtick.labelsize': 'medium',
                 'ytick.labelsize': 'medium',
                 'legend.fontsize': 'medium',
                 'figure.dpi': 160,
                 'savefig.dpi': 160,
                 'font.family': 'sans-serif',
                 'font.serif': ['Computer Modern Roman'],
                 # 'text.usetex': True,
                 'figure.figsize': (columnwidth_in, columnwidth_in / 2)}
