"""
This script will run all jupyter notebooks in order to test for errors.
"""
import sys
import os

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

os.chdir(os.path.dirname(sys.argv[0]))

notebooks = ('grids-and-coefficients.ipynb',
             'localized-spectral-analysis.ipynb',
             'gravity-and-magnetic-fields.ipynb',
             'plotting-maps.ipynb',
             'low-level-spherical-harmonic-analyses.ipynb',
             'advanced-localized-spectral-analysis.ipynb',
             'advanced-shcoeffs-and-shgrid-usage.ipynb',
             'spherical-harmonic-normalizations.ipynb',
             'advanced-shwindow-usage.ipynb',
             '3d-plots.ipynb')

if sys.version_info.major == 3:
    kname = 'python3'
else:
    raise ('Python version {:d} not supported.'.format(sys.version_info.major))

for i in range(len(notebooks)):
    with open(notebooks[i]) as f:
        nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=240, kernel_name=kname)

        print('Processing file {:s}'.format(notebooks[i]))
        ep.preprocess(nb, {'metadata': {'path': '.'}})
