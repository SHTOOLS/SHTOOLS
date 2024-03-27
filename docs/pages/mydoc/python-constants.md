---
title: "Constants"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-constants.html
summary: pyshtools provides easy access to many research-grade planetary constants.
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

## Constants

The *constants* subpackage defines physical constants related to the terrestrial planets and moons. Each constant is an instance of an [astropy](http://docs.astropy.org/en/stable/constants/index.html) `Constant` class, which has the attributes `name`, `value`, `uncertainty`, `unit`, and `reference`.

Each body can have several attributes, including
* `gm`,
* `mass`,
* `mean_radius` (aliased as `r`),
* `volume_equivalent_radius`,
* `volume`
* `mean_density`
* `gravity_mean_radius`
* `angular_velocity`
* `rotational_period`
* `orbit_angular_velocity`
* `orbit_semimajor_axis`
* `orbit_eccentricity`
* `orbit_inclination`
* `orbit_period`

Additional parameters are defined when appropriate. To see all information about an individual constant, it is only necessary to use the print function:
```python
In [1]: print(pysh.constants.Mars.mean_radius)
  Name   = Mean radius of Mars
  Value  = 3389500.0
  Uncertainty  = 0.0
  Unit  = m
  Reference = MOLA_shape: Wieczorek, M. (2024). Spherical harmonic models of the shape of Mars (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10794059
```
To use the value of a constant in a calculation, such as in this simple calculation of the circumference in kilometers, it is only necessary to access its `value` attribute:
```python
In [2]: 2 * np.pi * pysh.constants.Mars.mean_radius.value / 1000
21296.856598685208
```
To convert to a different set of units, one only needs to use the `to` method. For example, this converts the rotational period of Callisto from seconds to days
```python
In [3]: pysh.constants.Callisto.rotational_period.to('day')
<Quantity 16.68901797 d>
```
and this computes the gravitational acceleration on the mean planetary radius of
Mercury and then returns the value as a simple float in mGals:
```python
In [4]: (Mercury.gm / Mercury.mean_radius**2).to_value('mGal')
370218.70697392424
```
Physical constants from the *Committee on Data for Science and Technology* are provided in the submodule `codata`, and a few of these (such as `G` and `mu0`) are referenced in the main constants namespace.
