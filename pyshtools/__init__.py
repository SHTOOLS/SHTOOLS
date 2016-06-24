"""
Python Wrapper for the SHTOOLS library by Mark Wieczorek

see   http://shtools.ipgp.fr/
or    https://github.com/SHTOOLS/SHTOOLS

The python wrapper was written by: Matthias Meschede, Mark Wieczorek, 2014
"""

from __future__ import absolute_import, division, print_function

__version__ = '3.1'

#---- some miniature python functions: ----
def PlmIndex(l,m):
    return (l * (l + 1)) // 2 + m

def YilmIndexVector(l,m):
    return l**2 + (i-1)*l + m

#---- load docstring functions ----
def load_documentation():
    """
    Fills the modules __doc__ strings with a useful documentation that was
    generated at compile time
    """

    import os
    from . import _SHTOOLS

    # bind python functions to SHTOOLS
    _SHTOOLS.PlmIndex = PlmIndex
    _SHTOOLS.YilmIndexVector = YilmIndexVector

    print('Loading SHTOOLS documentation')
    pydocfolder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'doc'))
    for name, func in _SHTOOLS.__dict__.items():
        if callable(func):
            try:
                path = os.path.join(pydocfolder, name.lower() + '.doc')

                with open(path) as pydocfile:
                    pydoc = pydocfile.read()

                func.__doc__ = pydoc
            except IOError as msg:
                print(msg)

# load documentation that was generated at compile time (.doc files)
load_documentation()

# import planetary constants into module namespace
from . import _constant
constant = _constant.planetsconstants

# import all functions into module namespace
from ._SHTOOLS import *

# import class interface into namespace
from .classes import SHCoeffs, SHGrid, SHWindow

