"""
pyshtools constants.

To view information about a pyshtools constant, use

    pyshtools.constant.constantname.info()
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import os as _os
import numpy as _np

# ---- Import SHTOOLS constants into _constant
from .. import _constant


# -----------------------------------------------------------------------------
# --- Define a subclass of numpy.ndarray that adds an info() method for
# --- displaying documentation about a pyshtools constant.
# -----------------------------------------------------------------------------
class ndarrayinfo(_np.ndarray):
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


for _name, _value in _constant.planetsconstants.__dict__.items():
    locals()[_name] = _value.view(ndarrayinfo)


# -----------------------------------------------------------------------------
# ---- Fill the pyshtools constant infostrings with documentation from
# ---- external files. The doc files are generated during intitial compilation
# ---- of pyshtools from md formatted text files.
# -----------------------------------------------------------------------------
_pydocfolder = _os.path.abspath(_os.path.join(
                   _os.path.split(_os.path.dirname(__file__))[0], 'doc'))

for _name in _constant.planetsconstants.__dict__.keys():
    try:
        _path = _os.path.join(_pydocfolder, 'constant_' + _name.lower() +
                              '.doc')

        with open(_path, 'rb') as _pydocfile:
            _pydoc = _pydocfile.read().decode('utf-8')

        setattr(locals()[_name], '_infostring', _pydoc)

    except IOError as msg:
        print(msg)
