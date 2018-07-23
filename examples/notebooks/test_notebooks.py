"""
This script will run all jupyter notebooks in order to test for errors.
"""
from __future__ import print_function

import sys
import os

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

os.chdir(os.path.dirname(sys.argv[0]))

notebooks = ('Introduction-1.ipynb',
             'Introduction-2.ipynb',
             'Introduction-3.ipynb',
             'tutorial_1.ipynb',
             'tutorial_2.ipynb',
             'tutorial_3.ipynb',
             'tutorial_4.ipynb',
             'tutorial_5.ipynb',
             'tutorial_6.ipynb')

if sys.version_info.major == 2:
    kname = 'python2'
elif sys.version_info.major == 3:
    kname = 'python3'
else:
    raise ('Python version {:d} not supported.'.format(sys.version_info.major))

print('Python kernel name = {:s}'.format(kname))

for i in range(len(notebooks)):
    with open(notebooks[i]) as f:
        nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=240, kernel_name=kname)

        print('Processing file {:s}'.format(notebooks[i]))
        ep.preprocess(nb, {'metadata': {'path': '.'}})
