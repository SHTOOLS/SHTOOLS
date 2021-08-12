# defaults
_preferred_backend = "shtools"
from . import shtools as _preferred_backend_module


def preferred_backend():
    return _preferred_backend


def preferred_backend_module():
    return _preferred_backend_module


def backend_module(backend=None, nthreads=None):
    """Return the backend module."""
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
