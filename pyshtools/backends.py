"""
backends

This module defines functions for accessing and setting optional parameters of
the backends used for the spherical harmonic transforms in pyshtools.

Supported backends
------------------
shtools (default)            Spherical Harmonic Tools (Fortran 95), installed
                             as part of pyshtools.
ducc0                        Distinctly Useful Code Collection (C++17).

Functions
---------
select_preferred_backend()   Set the preferred backend module and options.
preferred_backend()          Return the name of the current preferred backend.
preferred_backend_module()   Return a reference to the preferred backend
                             module.
backend_module()             Return a reference to the specified backend
                             module.
"""

# defaults
from . import shtools as _preferred_backend_module
_preferred_backend = "shtools"


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
    return _preferred_backend_module


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
    backend = backend.lower()
    if backend == "shtools":
        from . import shtools

        return shtools
    elif backend == "ducc":
        from .wrappers import ducc0_wrapper
        if not ducc0_wrapper.available():
            raise ImportError('"ducc" backend requested, but not installed.')
        if nthreads is not None:
            ducc0_wrapper.set_nthreads(nthreads)
        return ducc0_wrapper
    elif backend is None:
        return preferred_backend_module()
    else:
        print("Unknown backend '{}' requested.".format(backend))
        raise RuntimeError


def select_preferred_backend(backend="shtools", nthreads=None):
    """
    Select the preferred backend module used for the spherical harmonic
    transforms in pyshtools.

    Usage
    -----
    select_preferred_backend_module([backend, nthreads])

    Parameters
    ----------
    backend : str, optional, default = 'shtools'
        Name of the preferred backend, either 'shtools' or 'ducc'.
    nthreads : int, optional, default = 1
        Number of threads to use for the 'ducc' backend. Setting this parameter
        to 0 will use as many threads as there are hardware threads on the
        system.
    """
    global _preferred_backend, _preferred_backend_module
    backend = backend.lower()
    if backend == "shtools":
        _preferred_backend = backend
        from . import shtools

        _preferred_backend_module = shtools
    elif backend == "ducc":
        try:
            import ducc0

            major, minor, patch = ducc0.__version__.split(".")
            if int(major) < 1 and int(minor) < 15:
                print(
                    "ducc0 installation found, but it is too old. "
                    "Need at least version 0.15"
                )
                raise RuntimeError
        except:
            print(
                "DUCC backend requested, but the relevant package cannot be "
                "imported. Leaving backend unchanged."
            )
        _preferred_backend = backend
        from .wrappers import ducc0_wrapper

        _preferred_backend_module = ducc0_wrapper
        if nthreads is not None:
            ducc0_wrapper.set_nthreads(nthreads)
    else:
        print("Unknown backend '{}' requested.".format(backend))
        raise RuntimeError


__all__ = ['preferred_backend', 'preferred_backend_module', 'backend_module',
           'select_preferred_backend']
