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
from . import ducc0_wrapper
from . import shtools

# defaults
_preferred_backend = "shtools"
_available_backends = {"shtools": shtools}

if ducc0_wrapper.available():
    _preferred_backend = "ducc"
    _available_backends["ducc"] = ducc0_wrapper


def preferred_backend():
    """
    Return the name of the current preferred backend used for the spherical
    harmonic transforms in pyshtools as a string.

    Usage
    -----
    preferred_backend()
    """
    return _preferred_backend


def preferred_backend_module():
    """
    Return a reference to the current preferred backend module used for the
    spherical harmonic transforms in pyshtools.

    Usage
    -----
    module = preferred_backend_module()
    """
    return _available_backends[_preferred_backend]


def backend_module(backend=None, nthreads=None):
    """
    Return the specified backend module used for the spherical harmonic
    transforms in pyshtools.

    Usage
    -----
    module = backend_module([backend, nthreads])

    Parameters
    ----------
    backend : str, optional, default = preferred_backend()
        Name of the preferred backend, either 'shtools' or 'ducc'.
    nthreads : int, optional, default = 1
        Number of threads to use for the 'ducc' backend. Setting this parameter
        to 0 will use as many threads as there are hardware threads on the
        system.
    """
    if backend is None:
        return preferred_backend_module()
    backend = backend.lower()
    if backend not in _available_backends:
        print("Requested backend '{}' not available.".format(backend))
        raise RuntimeError
    if backend == "ducc" and nthreads is not None:
        ducc0_wrapper.set_nthreads(nthreads)
    return _available_backends[backend]


def select_preferred_backend(backend=_preferred_backend, nthreads=None):
    """
    Select the preferred backend module used for the spherical harmonic
    transforms in pyshtools.

    Usage
    -----
    select_preferred_backend([backend, nthreads])

    Parameters
    ----------
    backend : str, optional, default = 'ducc'
        Name of the preferred backend, either 'shtools' or 'ducc'.
    nthreads : int, optional, default = 1
        Number of threads to use for the 'ducc' backend. Setting this parameter
        to 0 will use as many threads as there are hardware threads on the
        system.
    """
    global _preferred_backend
    backend = backend.lower()
    if backend in _available_backends:
        _preferred_backend = backend
        if backend == "ducc" and nthreads is not None:
            ducc0_wrapper.set_nthreads(nthreads)
        # inject functions
        from ..expand import inject_backend_specific_functions_for_expand
        inject_backend_specific_functions_for_expand()
        from ..rotate import inject_backend_specific_functions_for_rotate
        inject_backend_specific_functions_for_rotate()
    else:
        print("Requested backend '{}' not available.".format(backend))
        raise RuntimeError


__all__ = ['preferred_backend', 'preferred_backend_module', 'backend_module',
           'select_preferred_backend']
