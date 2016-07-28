"""
pyshtools
=========

pyshtools is an archive of scientific routines that can be used to
perform spherical harmonic transforms and reconstructions, rotations
of data expressed in spherical harmonics, and multitaper spectral
analyses on the sphere.

This module makes use of Python-wrapped Fortran 95 routines. For
further information, consult the web documentation at

   http://shtools.ipgp.fr/

and the GitHub project page at

   https://github.com/SHTOOLS/SHTOOLS
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

__version__ = '3.3-beta'
__author__ = 'SHTOOLS developers'

import os as _os
import numpy as _np

# ---- Import wrapped SHTOOLS functions into _SHTOOLS
from . import _SHTOOLS

# ---- Import SHTOOLS constants into _constant
from . import _constant

# ---- Import classes into pyshtools namespace
from .shclasses import SHCoeffs, SHGrid, SHWindow


# ---------------------------------------------------------------------
# --- Define two Python functions that replace their Fortran
# --- equivalents that use different indexing conventions, then
# --- bind these function to _SHTOOLS.
# ---------------------------------------------------------------------
def PlmIndex(l, m):
    return (l * (l + 1)) // 2 + m


def YilmIndexVector(i, l, m):
    return l**2 + (i - 1) * l + m

_SHTOOLS.PlmIndex = PlmIndex
_SHTOOLS.YilmIndexVector = YilmIndexVector


# ---------------------------------------------------------------------
# --- Define a subclass of numpy.ndarray that adds an info() method for
# --- displaying documentation about a pyshtools constant. Then define
# --- ConstantClass that holds these objects.
# ---------------------------------------------------------------------
class _ndarrayinfo(_np.ndarray):
    """
    To view information about a pyshtools constant, use

    pyshtools.constant.constantname.info()
    """
    def __new__(cls, input_array, infostring=None):
        # Input array is an already formed ndarray instance
        obj = _np.asarray(input_array).view(cls)
        obj._infostring = infostring
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._infostring = getattr(obj, '_infostring', None)

    def info(self):
        """
        To view information about a pyshtools constant, use

        pyshtools.constant.constantname.info()
        """
        print(self._infostring, end='')


class _ConstantClass():
    """
    This class is filled with the pyshtools constants
    To view information about a pyshtools constant, use

    pyshtools.constant.constantname.info()
    """
    pass

constant = _ConstantClass()

for _name, _value in _constant.planetsconstants.__dict__.items():
    setattr(constant, _name, _value.view(_ndarrayinfo))


# ---------------------------------------------------------------------
# ---- Fill the pyshtools module doc strings and pyshtools constant
# ---- infostrings with documentation from external files. The doc files
# ---- are generated during intitial compilation of pyshtools from md
# ---- formatted text files.
# ---------------------------------------------------------------------
print('Loading SHTOOLS -- version', __version__)

_pydocfolder = _os.path.abspath(_os.path.join(_os.path.dirname(__file__),
                                              'doc'))

for _name, _func in _SHTOOLS.__dict__.items():
    if callable(_func):
        try:
            _path = _os.path.join(_pydocfolder, _name.lower() + '.doc')

            with open(_path) as _pydocfile:
                _pydoc = _pydocfile.read()

            _func.__doc__ = _pydoc
        except IOError as msg:
            print(msg)

for _name in _constant.planetsconstants.__dict__.keys():
    try:
        _path = _os.path.join(_pydocfolder, 'constant_' + _name.lower() +
                              '.doc')

        with open(_path) as _pydocfile:
            _pydoc = _pydocfile.read()

        setattr(getattr(constant, _name), '_infostring', _pydoc)

    except IOError as msg:
        print(msg)


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['_ndarrayinfo', '_ConstantClass', 'constant', 'shclasses']
__all__ += ['SHCoeffs', 'SHGrid', 'SHWindow']

# --- Import all functions into pyshtools namespace ----
for _name, _func in _SHTOOLS.__dict__.items():
    if callable(_func):
        __all__.append(_name)
        globals()[_name] = _func
