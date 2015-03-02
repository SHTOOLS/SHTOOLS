"""
Python Wrapper for the SHTOOLS library by Mark Wieczorek

Authors: Matthias Meschede, Mark Wieczorek, 2014
"""

def load_documentation():
    import os
    from . import _SHTOOLS
    print('loading shtools documentation')
    pydocfolder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'doc'))
    for name,func in _SHTOOLS.__dict__.items():
        if callable(func):
            try:
                path = os.path.join(pydocfolder,name.lower()+'.doc')

                pydocfile = open(path)
                pydoc = pydocfile.read()
                pydocfile.close()

                func.__doc__ = pydoc 
            except IOError as msg:
                print(msg)

#try to load documentation if interactive
import __main__ as main
if not hasattr(main, '__file__'):
    load_documentation()

#import into main namespace
from . import _constants
constants = _constants.planetsconstants

from _SHTOOLS import *
