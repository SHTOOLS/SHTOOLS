---
title: "Backends"
keywords: spherical harmonics software package, backends, shtools, ducc
sidebar: mydoc_sidebar
permalink: backends.html
summary: pyshtools supports multiple backends for performing the most computationaly intensive mathematical operations.
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
</style>


## pyshtools backends

*pyshtools* is a collection of utilities that facilitates working with data expressed in spherical harmonics. Though *pyshtools* is accessed through python, many of the mathematical operations that are called are computationally demanding and are instead perfomed by external libraries. In order to access the functions of these external libraries, an intermediary python wrapper function is used that converts the python call and parameters into the form expected by the library. In principle, any pre-existing numerical software package could be used as a numerical backend to *pyshtools* by the creation of such wrapper functions.

The following software packages are supported as numerical backends to *pysthools*:

| Class name | Description |
| ---------- | ----------- |
| [shtools](index-fortran.html) (default) | Spherical Harmonic Tools (Fortran 95, installed as part of pyshtools) |
| [ducc](https://gitlab.mpcdf.mpg.de/mtr/ducc) | Distinctly Useful Code Collection (C++17) |

The *shtools* package (written in Fortran 95) suppports all the functionality of *pyshtools*. Other packages only implement a subset of the *shtools* routines, and whenever an implementation of a function does not exist, the *shtools* function is used as a fallback.

There are two ways to control which backend is used. The first method is to set the preferred backend that will be used in all subsequent calls. As an example, here we select the *ducc* backend using the `select_preferred_backend()` function, and also set the number of threads to use to 4:
```python
In [1]: pysh.backends.select_preferred_backend(backend='ducc', nthreads=4)
```
The current preferred backend can be inspected by using the `preferred_backend()` function:
```python
In [2]: pysh.backends.preferred_backend()
Out[2]: 'ducc'
```

Alternatively, if one has an initialized *pyshtools* class instance (such as from `SHGrid` or `SHCoeffs`), the backend can be specified for an single function call by use of the optional `backend` parameter. In this example, `grid` is expanded into spherical harmonic coefficients using the same backend parameters as above:
```python
In [3]: clm = grid.expand(backend='ducc', nthreads=4)
```
Following such a call, the preferred backend remains unchanged.


## Supported backend key characteristics

### shtools

* Supports all *pyshtools* functionality, including localized spectral analyses, Slepian analyses, and operations on gravity and magnetic field data.
* Spherical harmonic transforms and reconstructions are accurate up to degree 2800.
* Spherical harmonic rotations are accurate to about degree 1200.

### ducc
* Spherical harmonic transforms and rotations are more than 10 times faster than the default *shtools* routines.
* Supports the use of multiple threads to speed up computations.
* Spherical harmonic transforms are accurate beyond degree 25,000.
* Does not implement functions involving localized spectral analyses, Slepian analyses, nor gravity and magnetic field data.
* The native routines make use orthonormalized spherical harmonics that exclude the Condon-Shortley phase factor: Other normalizations are supported, but require preprocessing in python.


## Speed comparisons

The speeds of the spherical harmonic transforms using the above backends were tested for both real and complex data using equidistant 'DH2' grids. The amount of time in seconds required to perform individually the forward and inverse operations is plotted in the figure below as a function of the spherical harmonic bandwidth of the function. These calculations were performed on a 2018 MacBook Pro with a 2.7 GHz Quad-Core Intel i7 processor and 16 GB of memory.

{% include image.html file="backend-timing-tests.png" alt="Timing tests" caption="Figure 1. Time to perform the reconstruction of a function from its spherical harmonic coefficients (solid lines) and the spherical harmonic transform of the gridded function (dashed lines). The spherical harmonic transforms made use of orthonormalized spherical harmonics that exclude the Condon-Shortley phase factor (which is the native implementation in the 'ducc' backend)." %}

For the 'shtools' backend, the transform time is on the order of one second for degrees close to 700 and just under a minute for degree 2800. The 'ducc' backend has transform speeds that are about an order of magnitude faster than the 'shtools' backend. Transforms are about 1 second for degrees near 2000, and between 1 and 3 minutes for degrees near 10,000. In general, the complex routines are slower than then real routines by a factor close to 2.
