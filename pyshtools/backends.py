_preferred_backend = "SHTOOLS"

# NOTE: it would be possible to pass the number of threads to use as an argument
# here. Currently I'm using the value stored in the OMP_NUM_THREADS environment
# variable, but this approach may be more flexible.
def select_ducc_backend(nthreads=None):
    global _preferred_backend
    try:
        import ducc0
        major, minor, patch = ducc0.__version__.split('.')
        if int(major) < 1 and int(minor) < 14:
            print("ducc0 installation found, but it is too old. "
                  "Need at least version 0.14")
            raise RuntimeError
        _preferred_backend = "DUCC"
        if nthreads is not None:
            from .wrappers import ducc0_wrapper
            ducc0_wrapper.set_nthreads(nthreads)
    except:
        print("DUCC backend requested, but the relevant package cannot be "
              "imported. Leaving backend unchanged.")

def select_shtools_backend():
    global _preferred_backend
    _preferred_backend = "SHTOOLS"


def preferred_backend():
    return _preferred_backend
