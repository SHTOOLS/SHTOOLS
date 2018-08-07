"""
A collection of matplotlib style sheets for creating publication quality
graphics.

Usage:

    import matplotlib.pyplot as plt
    import pyshtools.utils.figstyle as figstyle

    # use a single style file
    plt.style.use(figstyle.shtools)
    # combine multiple style files
    plt.style.use([figstyle.shtools, figstyle.half])
    # apply style temporarily
    with plt.style.context(figstyle.shtools):

Styles:

    shtools :        Core style file for pyshtools.
    full :           Set the figure width and height for figures that span the
                     entire width of a page.
    threequarters :  Set the figure width and height for figures that span
                     three-quarters of a page (default for shtools style).
    half :           Set the figure width and height for figures that span
                     half a page.
    map :            Set parameters for plotting global maps.
"""

width = 7.48031  # For Wiley-Blackwell AGU publications
aspect_ratio = 4 / 3

shtools = {
    'font.size': 10,
    'font.family': 'sans-serif',
    'font.sans-serif' : ['Myriad Pro', 'Arial', 'Helvetica'],
    'font.serif': ['Times'],
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'savefig.dpi': 600,
    'figure.dpi': 100,
    'figure.figsize' : (width * 3. / 4., width * 3. / 4. / aspect_ratio),
    'text.usetex' : False,
    'axes.linewidth' : 1,
    'grid.linewidth' : 0.3,
    'grid.color' : 'k',
    'grid.linestyle' : '--',
    'xtick.major.width' : 0.6,
    'ytick.major.width' : 0.6,
    'xtick.minor.width' : 0.6,
    'xtick.minor.width' : 0.6,
    'lines.linewidth' : 1.5,
    'legend.framealpha' : 1.,
    'legend.edgecolor' : 'k'
    }

full = {
    'figure.figsize' : (width, width / aspect_ratio)
    }

threequarters = {
    'figure.figsize' : (width * 3. / 4., width * 3. / 4. / aspect_ratio)
    }

half = {
    'figure.figsize' : (width / 2., width / 2. / aspect_ratio)
    }

map = {
    'grid.linestyle' : '--',
    'grid.linewidth' : 0.25,
    'grid.color' : 'k'
    }


#font.sans-serif     : DejaVu Sans, Bitstream Vera Sans, Computer Modern Sans Serif, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica, Avant Garde, sans-serif


# 'lines.linewidth' : 1.5 def mpl, only for data lines
#font.family         : sans-serif
#font.style          : normal
#font.variant        : normal
#font.weight         : medium
#font.stretch        : normal
#text.color          : black

#axes.linewidth      : 0.8     # edge linewidth
#axes.grid           : False   # display grid or not

#axes.titlepad       : 6.0     # pad between axes and title in points
#axes.labelpad       : 4.0     # space between label and axis
#axes.labelweight    : normal  # weight of the x and y labels
#axes.labelcolor     : black

#axes.formatter.limits : -7, 7 # use scientific notation if log10
                               # of the axis range is smaller than the
                               # first or larger than the second
                               
#axes.xmargin        : .05  # x margin.  See `axes.Axes.margins`
#axes.ymargin        : .05  # y margin See `axes.Axes.margins`

#xtick.top            : False   # draw ticks on the top side
#xtick.bottom         : True   # draw ticks on the bottom side
#xtick.major.size     : 3.5      # major tick size in points
#xtick.minor.size     : 2      # minor tick size in points
#xtick.major.width    : 0.8    # major tick width in points
#xtick.minor.width    : 0.6    # minor tick width in points
#xtick.major.pad      : 3.5      # distance to major tick label in points
#xtick.minor.pad      : 3.4      # distance to the minor tick label in points
#xtick.color          : k      # color of the tick labels
#xtick.labelsize      : medium # fontsize of the tick labels
#xtick.direction      : out    # direction: in, out, or inout
#xtick.minor.visible  : False  # visibility of minor ticks on x-axis
#xtick.major.top      : True   # draw x axis top major ticks
#xtick.major.bottom   : True   # draw x axis bottom major ticks
#xtick.minor.top      : True   # draw x axis top minor ticks
#xtick.minor.bottom   : True   # draw x axis bottom minor ticks

#ytick.left           : True   # draw ticks on the left side
#ytick.right          : False  # draw ticks on the right side
#ytick.major.size     : 3.5      # major tick size in points
#ytick.minor.size     : 2      # minor tick size in points
#ytick.major.width    : 0.8    # major tick width in points
#ytick.minor.width    : 0.6    # minor tick width in points
#ytick.major.pad      : 3.5      # distance to major tick label in points
#ytick.minor.pad      : 3.4      # distance to the minor tick label in points
#ytick.color          : k      # color of the tick labels
#ytick.labelsize      : medium # fontsize of the tick labels
#ytick.direction      : out    # direction: in, out, or inout
#ytick.minor.visible  : False  # visibility of minor ticks on y-axis
#ytick.major.left     : True   # draw y axis left major ticks
#ytick.major.right    : True   # draw y axis right major ticks
#ytick.minor.left     : True   # draw y axis left minor ticks
#ytick.minor.right    : True   # draw y axis right minor ticks

#grid.color       :   b0b0b0    # grid color
#grid.linestyle   :   -         # solid
#grid.linewidth   :   0.8       # in points
#grid.alpha       :   1.0       # transparency, between 0.0 and 1.0

### Legend
#legend.loc           : best
#legend.frameon       : True     # if True, draw the legend on a background patch
#legend.framealpha    : 0.8      # legend patch transparency
#legend.facecolor     : inherit  # inherit from axes.facecolor; or color spec
#legend.edgecolor     : 0.8      # background patch boundary color
#legend.fancybox      : True     # if True, use a rounded box for the
                                 # legend background, else a rectangle
#legend.shadow        : False    # if True, give background a shadow effect
#legend.numpoints     : 1        # the number of marker points in the legend line
#legend.scatterpoints : 1        # number of scatter points
#legend.markerscale   : 1.0      # the relative size of legend markers vs. original
#legend.fontsize      : medium
# Dimensions as fraction of fontsize:
#legend.borderpad     : 0.4      # border whitespace
#legend.labelspacing  : 0.5      # the vertical space between the legend entries
#legend.handlelength  : 2.0      # the length of the legend lines
#legend.handleheight  : 0.7      # the height of the legend handle
#legend.handletextpad : 0.8      # the space between the legend line and legend text
#legend.borderaxespad : 0.5      # the border between the axes and legend edge
#legend.columnspacing : 2.0      # column separation

#figure.titlesize : large      # size of the figure title (Figure.suptitle())
#figure.titleweight : normal   # weight of the figure title
#figure.figsize   : 6.4, 4.8   # figure size in inches
#figure.dpi       : 100      # figure dots per inch
#figure.facecolor : white   # figure facecolor; 0.75 is scalar gray
#figure.edgecolor : white   # figure edgecolor
#figure.autolayout : False  # When True, automatically adjust subplot
                            # parameters to make the plot fit the figure

#figure.subplot.left    : 0.125  # the left side of the subplots of the figure
#figure.subplot.right   : 0.9    # the right side of the subplots of the figure
#figure.subplot.bottom  : 0.11    # the bottom of the subplots of the figure
#figure.subplot.top     : 0.88    # the top of the subplots of the figure
#figure.subplot.wspace  : 0.2    # the amount of width reserved for blank space between subplots,
                                 # expressed as a fraction of the average axis width
#figure.subplot.hspace  : 0.2    # the amount of height reserved for white space between subplots,
                                 # expressed as a fraction of the average axis height

#image.aspect : equal             # equal | auto | a number
#image.interpolation  : nearest   # see help(imshow) for options
#image.cmap   : viridis           # A colormap name, gray etc...
#image.lut    : 256               # the size of the colormap lookup table
#image.origin : upper             # lower | upper
#image.resample  : True
#image.composite_image : True     # When True, all the images on a set of axes are
                                  # combined into a single composite image before
                                  # saving a figure as a vector graphics file,
                                  # such as a PDF.

#savefig.dpi         : figure   # figure dots per inch or 'figure'
#savefig.facecolor   : white    # figure facecolor when saving
#savefig.edgecolor   : white    # figure edgecolor when saving
#savefig.format      : png      # png, ps, pdf, svg
#savefig.bbox        : standard # 'tight' or 'standard'.
                                # 'tight' is incompatible with pipe-based animation
                                # backends but will workd with temporary file based ones:
                                # e.g. setting animation.writer to ffmpeg will not work,
                                # use ffmpeg_file instead
#savefig.pad_inches  : 0.1      # Padding to be used when bbox is set to 'tight'
#savefig.jpeg_quality: 95       # when a jpeg is saved, the default quality parameter.
#savefig.directory   : ~        # default directory in savefig dialog box,
                                # leave empty to always use current working directory
#savefig.transparent : False    # setting that controls whether figures are saved with a
                                # transparent background by default

#pdf.compression   : 6 # integer from 0 to 9
                       # 0 disables compression (good for debugging)

#figure.autolayout : False     ## When True, automatically adjust subplot
                               ## parameters to make the plot fit the figure
                               ## using `tight_layout`
                               
