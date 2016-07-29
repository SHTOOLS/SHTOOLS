"""
pyshtools constants.

To view information about a pyshtools constant, use

    pyshtools.constant.constantname.info()
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as np

# ---- Import SHTOOLS constants into _constant
from . import _constant


# ---------------------------------------------------------------------
# --- Define a subclass of numpy.ndarray that adds an info() method for
# --- displaying documentation about a pyshtools constant.
# ---------------------------------------------------------------------
class ndarrayinfo(np.ndarray):
    """
    To view information about a pyshtools constant, use

    pyshtools.constant.constantname.info()
    """
    def __new__(cls, input_array, infostring=None):
        # Input array is an already formed ndarray instance
        obj = np.asarray(input_array).view(cls)
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
    globals()[_name] = _value.view(ndarrayinfo)
