"""
pyshtools Spherical Harmonic Rotation Routines.

This subpackage of pyshtools defines the following functions:

djpi2             Compute the rotation matrix d(pi/2) used in rotating data
                  expressed in spherical harmonics.
SHRotateCoef      Determine the spherical harmonic coefficients of a complex
                  function rotated by three Euler angles.
SHRotateRealCoef  Determine the spherical harmonic coefficients of a real
                  function rotated by three Euler angles.
"""
from ..backends.shtools import djpi2
from ..backends.shtools import SHRotateCoef
from ..backends import backend_module, select_preferred_backend


def inject_backend_specific_functions_for_rotate():
    mod = backend_module()
    global SHRotateRealCoef
    SHRotateRealCoef = mod.SHRotateRealCoef


# trigger the injection of the backend-specific functions
select_preferred_backend()

# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['djpi2', 'SHRotateCoef', 'SHRotateRealCoef']
