---
title: "SHGravGrid class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shgravgrid.html
summary: A class for global gridded gravitational field data.
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
| `x = SHGravCoeffs.expand()` | Initialize using an SHGravCoeffs class instance. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `rad` | SHGrid class instance of the radial component of the gravitational acceleration evaluated on an ellipsoid. |
| `theta` | SHGrid class instance of the theta component of the gravitational acceleration evaluated on an ellipsoid. |
| `phi` | SHGrid class instance of the phi component of the gravitational acceleration evaluated on an ellipsoid. |
| `total` | SHGrid class instance of the total gravitational acceleration with the normal gravity removed on an ellipsoid. |
| `pot` | SHGrid class instance of the gravitational potential evaluated on an ellipsoid. |
| `gm` | Gravitational constant time the mass of the body. |
| `a` | Semimajor axis of the reference ellipsoid. |
| `f` | Flattening of the reference ellipsoid, f = (a - b) / a. |
| `omega` | Angular rotation rate of the body. |
| `normal_gravity` | True if the normal gravity is removed from the total gravitational acceleration. |
| `lmax` | The maximum spherical harmonic degree resolvable by the grids. |
| `lmax_calc` | The maximum spherical harmonic degree of the gravitational potential used in creating the grids. |
| `units` | The units of the gridded gravity data. |
| `pot_units` | The units of the gridded gravitational potential data. |
| `nlat`, `nlon` | The number of latitude and longitude bands in the grids. |
| `sampling` | The longitudinal sampling scheme of the grids: either 1 for `nlon` = `nlat` or 2 for `nlon` = 2 * `nlat`. |
| `epoch` | The epoch time of the gravity model. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `plot()` | Plot all three components of the gravity field and the total gravity disturbance.|
| `plot_rad()` | Plot the radial component of the gravity field. |
| `plot_theta()` | Plot the theta component of the gravity field. |
| `plot_phi()` | Plot the phi component of the gravity field. |
| `plot_total()` | Plot the total gravity disturbance. |
| `plot_pot()` | Plot the gravitational potential. |
| `to_xarray()` | Return the gravity gridded data as an xarray DataSet. |
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the SHGravGrid instance. |
