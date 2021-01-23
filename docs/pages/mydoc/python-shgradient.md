---
title: "SHGradient class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shgradient.html
summary: A class for the horizontal components of the gradient of a scalar.
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:70%;
}
</style>

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = SHGradient.gradient()` | Initialize using an SHCoeffs class instance. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `theta` | SHGrid class instance of the theta component of the horizontal gradient. |
| `phi` |  SHGrid class instance of the phi component of the horizontal gradient. |
| `lmax` | The maximum spherical harmonic degree resolvable by the grids. |
| `lmax_calc` | The maximum spherical harmonic degree used in creating the grids. |
| `units` | The units of the gridded gradients. |
| `nlat`, `nlon` | The number of latitude and longitude bands in the grids. |
| `sampling` | The longitudinal sampling scheme of the grids: either 1 for `nlon` = `nlat` or 2 for `nlon` = 2 * `nlat`. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `plot()` | Plot the two components of the horizontal gradient. |
| `plot_theta()` | Plot the theta component of the horizontal gradient. |
| `plot_phi()` | Plot the phi component of the horizontal gradient. |
| `to_xarray()` | Return the gridded gradient data as an xarray DataSet. |
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the SHGradient instance. |
