"""
backends

This module defines functions for accessing and setting optional parameters of
the backends used for the spherical harmonic transforms in pyshtools.

Supported backends
------------------
shtools (fallback)           Spherical Harmonic Tools (Fortran 95), installed
                             as part of pyshtools.
ducc0 (default)              Distinctly Useful Code Collection (C++17).

Functions
---------
select_preferred_backend()   Set the preferred backend module and options.
preferred_backend()          Return the name of the current preferred backend.
preferred_backend_module()   Return a reference to the preferred backend
                             module.
backend_module()             Return a reference to the specified backend
                             module.
"""
from .backends import select_preferred_backend
from .backends import preferred_backend
from .backends import preferred_backend_module
from .backends import backend_module
from . import ducc0_wrapper  # noqa: F401
from . import shtools  # noqa: F401

del backends  # noqa: F821


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['select_preferred_backend', 'preferred_backend',
           'preferred_backend_module', 'backend_module']
