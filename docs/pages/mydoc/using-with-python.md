---
title: "Using pyshtools"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: using-with-python.html
summary: To use pyshtools in Python, it is only necessary to import the pyshtools module.
toc: true
folder: mydoc
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

To use pyshtools in Python, it is only necessary to execute this statement in the Python environment
```python
import pyshtools as pysh
```
This will load the following classes and subpackages into the `pysh` namespace:

| Class | Description |
| ----- | ----------- |
| [SHCoeffs](python-shcoeffs.html) | Class for spherical harmonic coefficients |
| [SHGrid](python-shgrid.html) | Class for global grids |
| [SHWindow](python-shwindow.html) | Class for localization windows |
| [Slepian](python-slepian.html) | Class for Slepian functions |
| [SHGravCoeffs](python-shgravcoeffs.html) | Class for gravitational potential spherical harmonic coefficients |
| [SHMagCoeffs](python-shmagcoeffs.html) | Class for magnetic potential spherical harmonic coefficients |


| Subpackage | Description |
| ---------- | ----------- |
| [shclasses](python-classes.html) | All classes and subclasses |
| [constants](python-datasets-constants.html) | Planetary constants |
| [datasets](python-datasets-constants.html) | Gravity, topography and magnetic field datasets |
| [legendre](python-legendre-functions.html) | Legendre functions |
| [expand](python-spherical-harmonic-transforms.html) | Spherical harmonic expansion routines |
| [shio](python-io.html) | Spherical harmonic I/O, storage, and conversion routines |
| [spectralanalysis](python-spectral-analysis.html) | Global and localized spectral analysis routines |
| [rotate](python-spherical-harmonic-rotations.html) | Spherical harmonic rotation routines |
| [gravmag](python-gravity-magnetics.html) | Gravity and magnetics routines |
| [utils](python-utilities.html) | Utilities |
| [backends](python-backends.html) | Functions to control which backend to use for the spherical harmonic transforms. |
| shtools | All Python wrapped Fortran 95 routines |


If you are using [iPython](https://ipython.org), which adds improved functionality to Python, the available pyshtools routines can be explored by typing
```python
pysh.[tab]
```
where `[tab]` is the tab key.

To read the documentation of a routine in iPython, such as the `expand()` method of the `SHCoeffs` class, enter
```python
pysh.SHCoeffs.expand?
```

To use the pyshtools matplotlib style parameters for publication quality graphics, input
```python
pysh.utils.figstyle()
```
This function takes optional parameters for specifying the screen resolution, the relative width of the figure, the physical width of the journal page in inches, and the aspect ratio of the figure.

