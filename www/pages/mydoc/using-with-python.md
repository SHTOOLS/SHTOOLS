---
title: "Using SHTOOLS in Python"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: using-with-python.html
summary: To use SHTOOLS in Python, it is only necessary to import the pyshtools module.
toc: true
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>


## The pyshtools module

To use SHTOOLS in Python, it is only necessary to execute this statement in the Python environment
```python
import pyshtools
```
This will load the following three classes and subpackages into the `pyshtools` namespace:

| Class | Description |
| ----- | ----------- |
| [SHCoeffs](pyshcoeffs.html) | A high level class for spherical harmonic coefficients |
| [SHGrid](pyshgrid.html) | A high level class for global grids |
| [SHWindow](pyshwindow.html) | A high level class for localization windows |


| Subpackage | Description |
| ---------- | ----------- |
| [shclasses](python-classes.html) | All pyshtools classes and subclasses |
| shtools | All Python wrapped Fortran 95 routines |
| [legendre](python-legendre-functions.html) | Legendre functions |
| [expand](python-spherical-harmonic-transforms.html) | Spherical harmonic expansion routines |
| [shio](python-io.html) | Spherical harmonic I/O, storage, and conversion routines |
| [spectralanalysis](python-spectral-analysis.html) | Global and localized spectral analysis routines |
| [rotate](python-spherical-harmonic-rotations.html) | Spherical harmonic rotation routines |
| [gravmag](python-gravity-magnetics.html) | Gravity and magnetics routines |
| constant | pyshtools constants |
| [utils](python-utilities.html) | Utilities |

If you are using [iPython](http://ipython.org), which adds improved functionality to Python, the available `pyshtools` routines can be explored by typing
```python
pyshtools.[tab]
```
where `[tab]` is the tab key.

## Documentation

To read the documentation of a routine in iPython, such as `MakeGridDH`, enter
```python
pyshtools.expand.MakeGridDH?
```
To read the info string of an SHTOOLS constant, such as `a_mars`, enter
```python
pyshtools.constant.a_mars.info()
```
Documentation for the Python functions used in SHTOOLS can also be accessed by their unix man pages, appending `py` to the name and using all lower case letters. As an example, to access the python `MakeGridDH` man page, use
```bash
man pymakegriddh
```
Alternatively, the man pages can be accessed from the *Python components* menu item on this web site.
