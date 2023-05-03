"""
pyshtools subpackage that includes all Python wrapped Fortran routines.
"""
import os as _os
import numpy as _np

# Import all wrapped SHTOOLS functions

# legendre
from .._SHTOOLS import PlmBar
from .._SHTOOLS import PlmBar_d1
from .._SHTOOLS import PlBar
from .._SHTOOLS import PlBar_d1
from .._SHTOOLS import PlmON
from .._SHTOOLS import PlmON_d1
from .._SHTOOLS import PlON
from .._SHTOOLS import PlON_d1
from .._SHTOOLS import PlmSchmidt
from .._SHTOOLS import PlmSchmidt_d1
from .._SHTOOLS import PlSchmidt
from .._SHTOOLS import PlSchmidt_d1
from .._SHTOOLS import PLegendreA
from .._SHTOOLS import PLegendreA_d1
from .._SHTOOLS import PLegendre
from .._SHTOOLS import PLegendre_d1

# expand
from .._SHTOOLS import SHExpandDH
from .._SHTOOLS import MakeGridDH
from .._SHTOOLS import SHExpandDHC
from .._SHTOOLS import MakeGridDHC
from .._SHTOOLS import SHGLQ
from .._SHTOOLS import SHExpandGLQ
from .._SHTOOLS import MakeGridGLQ
from .._SHTOOLS import SHExpandGLQC
from .._SHTOOLS import MakeGridGLQC
from .._SHTOOLS import GLQGridCoord
from .._SHTOOLS import SHExpandLSQ
from .._SHTOOLS import SHExpandWLSQ
from .._SHTOOLS import MakeGrid2D
from .._SHTOOLS import MakeGridPoint
from .._SHTOOLS import MakeGridPointC
from .._SHTOOLS import SHMultiply
from .._SHTOOLS import MakeGradientDH

# shio
from .._SHTOOLS import SHRead2
from .._SHTOOLS import SHRead2Error
from .._SHTOOLS import SHReadJPL
from .._SHTOOLS import SHReadJPLError
from .._SHTOOLS import SHCilmToCindex
from .._SHTOOLS import SHCindexToCilm
from .._SHTOOLS import SHCilmToVector
from .._SHTOOLS import SHVectorToCilm
from .._SHTOOLS import SHrtoc
from .._SHTOOLS import SHctor

# spectralanalysis
from .._SHTOOLS import SHAdmitCorr
from .._SHTOOLS import SHConfidence
from .._SHTOOLS import SHMultiTaperSE
from .._SHTOOLS import SHMultiTaperCSE
from .._SHTOOLS import SHLocalizedAdmitCorr
from .._SHTOOLS import SHReturnTapers
from .._SHTOOLS import SHReturnTapersM
from .._SHTOOLS import ComputeDm
from .._SHTOOLS import ComputeDG82
from .._SHTOOLS import SHFindLWin
from .._SHTOOLS import SHBiasK
from .._SHTOOLS import SHMTCouplingMatrix
from .._SHTOOLS import SHBiasAdmitCorr
from .._SHTOOLS import SHMTDebias
from .._SHTOOLS import SHMTVarOpt
from .._SHTOOLS import SHMTVar
from .._SHTOOLS import SHSjkPG
from .._SHTOOLS import SHMultiTaperMaskSE
from .._SHTOOLS import SHMultiTaperMaskCSE
from .._SHTOOLS import SHReturnTapersMap
from .._SHTOOLS import SHBiasKMask
from .._SHTOOLS import ComputeDMap
from .._SHTOOLS import Curve2Mask
from .._SHTOOLS import SHBias
from .._SHTOOLS import SphericalCapCoef
from .._SHTOOLS import SHRotateTapers
from .._SHTOOLS import SlepianCoeffs
from .._SHTOOLS import SlepianCoeffsToSH
from .._SHTOOLS import SHSCouplingMatrix
from .._SHTOOLS import SHSlepianVar
from .._SHTOOLS import SHSCouplingMatrixCap

# rotate
from .._SHTOOLS import djpi2
from .._SHTOOLS import SHRotateCoef
from .._SHTOOLS import SHRotateRealCoef

# gravmag
from .._SHTOOLS import MakeGravGridDH
from .._SHTOOLS import MakeGravGridPoint
from .._SHTOOLS import MakeGravGradGridDH
from .._SHTOOLS import MakeGeoidGridDH
from .._SHTOOLS import CilmPlusDH
from .._SHTOOLS import CilmMinusDH
from .._SHTOOLS import CilmPlusRhoHDH
from .._SHTOOLS import CilmMinusRhoHDH
from .._SHTOOLS import BAtoHilmDH
from .._SHTOOLS import BAtoHilmRhoHDH
from .._SHTOOLS import DownContFilterMA
from .._SHTOOLS import DownContFilterMC
from .._SHTOOLS import NormalGravity
from .._SHTOOLS import MakeMagGridDH
from .._SHTOOLS import MakeMagGridPoint
from .._SHTOOLS import MakeMagGradGridDH

# utils
from .._SHTOOLS import MakeCircleCoord
from .._SHTOOLS import MakeEllipseCoord
from .._SHTOOLS import Wigner3j
from .._SHTOOLS import DHaj

__all__ = ['PlmBar', 'PlmBar_d1', 'PlBar', 'PlBar_d1', 'PlmON', 'PlmON_d1',
           'PlON', 'PlON_d1', 'PlmSchmidt', 'PlmSchmidt_d1', 'PlSchmidt',
           'PlSchmidt_d1', 'PLegendreA', 'PLegendreA_d1', 'PLegendre',
           'PLegendre_d1', 'SHExpandDH', 'MakeGridDH', 'SHExpandDHC',
           'MakeGridDHC', 'SHGLQ', 'SHExpandGLQ', 'MakeGridGLQ',
           'SHExpandGLQC', 'MakeGridGLQC', 'GLQGridCoord', 'SHExpandLSQ',
           'SHExpandWLSQ', 'MakeGrid2D', 'MakeGridPoint', 'MakeGridPointC',
           'SHMultiply', 'SHRead2', 'SHRead2Error', 'SHReadJPL',
           'SHReadJPLError', 'SHCilmToVector', 'SHVectorToCilm',
           'SHCilmToCindex', 'SHCindexToCilm', 'SHrtoc', 'SHctor',
           'SHAdmitCorr', 'SHConfidence', 'SHMultiTaperSE', 'SHMultiTaperCSE',
           'SHLocalizedAdmitCorr', 'SHReturnTapers', 'SHReturnTapersM',
           'ComputeDm', 'ComputeDG82', 'SHFindLWin', 'SHBiasK',
           'SHMTCouplingMatrix', 'SHBiasAdmitCorr', 'SHMTDebias', 'SHMTVarOpt',
           'SHSjkPG', 'SHMultiTaperMaskSE', 'SHMultiTaperMaskCSE',
           'SHReturnTapersMap', 'SHBiasKMask', 'ComputeDMap', 'Curve2Mask',
           'SHBias', 'SphericalCapCoef', 'djpi2', 'SHRotateCoef',
           'SHRotateRealCoef', 'MakeGravGridDH', 'MakeGravGradGridDH',
           'MakeGeoidGridDH', 'CilmPlusDH', 'CilmMinusDH', 'CilmPlusRhoHDH',
           'CilmMinusRhoHDH', 'BAtoHilmDH', 'BAtoHilmRhoHDH',
           'DownContFilterMA', 'DownContFilterMC', 'NormalGravity',
           'MakeMagGridDH', 'MakeCircleCoord', 'MakeEllipseCoord', 'Wigner3j',
           'DHaj', 'MakeMagGradGridDH', 'SHRotateTapers', 'SlepianCoeffs',
           'SlepianCoeffsToSH', 'SHSCouplingMatrix', 'SHMTVar', 'SHSlepianVar',
           'SHSCouplingMatrixCap', 'MakeGravGridPoint', 'MakeMagGridPoint',
           'MakeGradientDH']

_fortran_functions = ['MakeGridPoint', 'MakeGridPointC', 'DownContFilterMA',
                      'DownContFilterMC', 'SHFindLWin', 'SHSjkPG',
                      'NormalGravity', 'SHConfidence', 'MakeGravGridPoint',
                      'MakeMagGridPoint']

_fortran_subroutines = list(set(__all__) - set(_fortran_functions))


# -----------------------------------------------------------------------------
#
#   Fill the module doc strings with documentation from external
#   files. The doc files are generated during intitial compilation of
#   pyshtools from md formatted text files.
#
# -----------------------------------------------------------------------------
_pydocfolder = _os.path.abspath(_os.path.join(
                   _os.path.split(_os.path.dirname(__file__))[0], 'doc'))

for _name in __all__:
    try:
        _path = _os.path.join(_pydocfolder, _name.lower() + '.doc')

        with open(_path, 'rb') as _pydocfile:
            _pydoc = _pydocfile.read().decode('utf-8')

        setattr(locals()[_name], '__doc__', _pydoc)

    except IOError as msg:
        print(msg)


# -----------------------------------------------------------------------------
#
#   Check the exit status of Fortran routines, raise exceptions, and
#   strip exitstatus from the Python return values.
#
# -----------------------------------------------------------------------------
class SHToolsError(Exception):
    pass


def _shtools_status_message(status):
    '''
    Determine error message to print when an SHTOOLS Fortran 95 routine exits
    improperly.
    '''
    if (status == 1):
        errmsg = 'Improper dimensions of input array.'
    elif (status == 2):
        errmsg = 'Improper bounds for input variable.'
    elif (status == 3):
        errmsg = 'Error allocating memory.'
    elif (status == 4):
        errmsg = 'File IO error.'
    else:
        errmsg = 'Unhandled Fortran 95 error.'
    return errmsg


def _raise_errors(func):
    def wrapped_func(*args, **kwargs):
        returned_values = func(*args, **kwargs)
        if returned_values[0] != 0:
            raise SHToolsError(_shtools_status_message(returned_values[0]))
        elif len(returned_values) == 2:
            return returned_values[1]
        else:
            return returned_values[1:]
    wrapped_func.__doc__ = func.__doc__
    return wrapped_func


for _func in _fortran_subroutines:
    locals()[_func] = _raise_errors(locals()[_func])


# -----------------------------------------------------------------------------
#
#   Vectorize some of the Fortran functions
#
# -----------------------------------------------------------------------------
MakeGridPoint = _np.vectorize(MakeGridPoint, excluded=[0])
MakeGridPointC = _np.vectorize(MakeGridPointC, excluded=[0])
DownContFilterMA = _np.vectorize(DownContFilterMA)
DownContFilterMC = _np.vectorize(DownContFilterMC)
NormalGravity = _np.vectorize(NormalGravity, excluded=[1, 2, 3, 4])
SHConfidence = _np.vectorize(SHConfidence)
