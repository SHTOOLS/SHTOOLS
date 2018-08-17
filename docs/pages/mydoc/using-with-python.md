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
| [SHCoeffs](python-shcoeffs.html) | Class for spherical harmonic coefficients |
| [SHGrid](python-shgrid.html) | Class for global grids |
| [SHWindow](python-shwindow.html) | Class for localization windows |
| [SHGravCoeffs](python-shgravcoeffs.html) | Class for gravitational potential spherical harmonic coefficients.|
| [SHGravGrid](python-shgravgrid.html) | Class for global gridded gravitational field data.|
| [SHGravTensor](python-shgravtensor.html) | Class for the gravity tensor and eigenvalues. |
| [SHGeoid](python-shgeoid.html) | Class for the geoid.|



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
| [constant](python-constants.html) | pyshtools constants |
| [utils](python-utilities.html) | Utilities |

To use the `pyshtools` matplotlib style parameters for publication quality graphics, input
```python
pyshtools.utils.figstyle()
```
This function takes optional parameters for specifying the screen resolution, the relative width of the figure, the physical width of the journal page in inches, and the aspect ratio of the figure.

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
Alternatively, the documentation can be accessed from the *Python components* menu item on this web site.

The `constant` subpackage defines physical constants related to the gravity, topography, and magnetic field of the terrestrial planets. Each of these is an instance of an [astropy](http://docs.astropy.org/en/stable/constants/index.html) `Constant` class, which has the attributes `name`, `value`, `uncertainty`, `unit`, and `reference`. To see all information about an individual constant, enter
```python
print(pyshtools.constant.r_mars)
```

