"""
pyshtools
=========

pyshtools is an archive of scientific routines that can be used to perform
spherical harmonic transforms and reconstructions, rotations of data expressed
in spherical harmonics, and multitaper spectral analyses on the sphere.

This module makes use of Python-wrapped Fortran 95 routines. For further
information, consult the web documentation at

   http://shtools.ipgp.fr/

and the GitHub project page at

   https://github.com/SHTOOLS/SHTOOLS
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function
import numpy as _np

# ---- Import class interface into pyshtools namespace ----
from .classes import SHCoeffs, SHGrid, SHWindow

from . import _SHTOOLS

__version__ = '3.3-beta'


# ---------------------------------------------------------------------
# --- Define two Python functions that replace their Fortran
# --- equivalents that use different indexing conventions.
# ---------------------------------------------------------------------
def PlmIndex(l, m):
    return (l * (l + 1)) // 2 + m


def YilmIndexVector(i, l, m):
    return l**2 + (i-1)*l + m


# ----------------------------------------------------------------------
# --- Define a subclass of ndarray that adds an info() method for
# --- displaying documentation about a pyshtools constant.
# ----------------------------------------------------------------------
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


def _load_documentation():
    """
    Fill the pyshtools module __doc__ strings and pyshtoools constant
    _infostrings with documentation from external files. The doc files
    are generated during intitial compilation of pyshtools from md
    formatted text files.
    """
    import os

    # import _SHTOOLS functions to add documentation
    # the _SHTOOLS functions are later imported into
    # the main namespace

    # bind Python functions to _SHTOOLS
    _SHTOOLS.PlmIndex = PlmIndex
    _SHTOOLS.YilmIndexVector = YilmIndexVector

    print('Loading SHTOOLS -- version', __version__)
    pydocfolder = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                  'doc'))
    for name, func in _SHTOOLS.__dict__.items():
        if callable(func):
            try:
                path = os.path.join(pydocfolder, name.lower() + '.doc')

                with open(path) as pydocfile:
                    pydoc = pydocfile.read()

                func.__doc__ = pydoc
            except IOError as msg:
                print(msg)

    for name in _constant.planetsconstants.__dict__.keys():
        try:
            path = os.path.join(pydocfolder, 'constant_' + name.lower() +
                                '.doc')

            with open(path) as pydocfile:
                pydoc = pydocfile.read()

            setattr(getattr(constant, name), '_infostring', pydoc)

        except IOError as msg:
            print(msg)


# --- Import all functions into pyshtools namespace ----
from ._SHTOOLS import *  # NOQA

# --- Import planetary constants into pyshtools namespace ----
from . import _constant  # NOQA

constant = _ConstantClass()

for _name, _value in _constant.planetsconstants.__dict__.items():
    constant._name = _value.view(_ndarrayinfo)

# ---- Load documentation that was generated at compile time ----
_load_documentation()

# ---- Define __all__ for use with 'from pyshtools import *' ----
__all__ = ['_ndarrayinfo', '_ConstantClass', 'constant', 'classes']
__all__ += ['SHCoeffs', 'SHRealCoefficients', 'SHComplexCoefficients',
            'SHGrid', 'DHGrid', 'GLQGrid', 'SHWindow', 'SHSymmetricWindow',
            'SHAsymmetricWindow']



for _name, _func in _SHTOOLS.__dict__.items():
    if callable(_func):
        __all__.append(_name)
