# backend_module()

Return the specified backend module used for the spherical harmonic
transforms in pyshtools.

# Usage

```python
module = backend_module([backend, nthreads])
```

# Parameters

**backend : str, optional, default = preferred_backend()**
:   Name of the preferred backend, either 'shtools' or 'ducc'.

**nthreads : int, optional, default = 1**
:   Number of threads to use for the 'ducc' backend. Setting this parameter
        to 0 will use as many threads as there are hardware threads on the
        system.
    