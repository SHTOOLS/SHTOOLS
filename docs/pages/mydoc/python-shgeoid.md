---
title: "SHGeoid class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shgeoid.html
summary: A class for gridded data of the geoid.
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
| `x = SHGravCoeffs.geoid()` | Initialize using an SHGravCoeffs class instance. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `geoid` | SHGrid class instance of the geoid. |
| `gm` | Gravitational constant time the mass of the body. |
| `potref` | Potential of the chosen geoid. |
| `a` | Semimajor axis of the reference ellipsoid. |
| `f` | Flattening of the reference ellipsoid, f = (a - b) / a. |
| `omega` | Angular rotation rate of the body. |
| `r` | Reference radius of the Taylor expansion. |
| `order` | Order of the Taylor expansion. |
| `units` | The units of the gridded data. |
| `lmax` | The maximum spherical harmonic degree resolvable by the grids. |
| `lmax_calc` | The maximum spherical harmonic degree of the gravitational potential used in creating the grids. |
| `nlat`, `nlon` | The number of latitude and longitude bands in the grids. |
| `sampling` | The longitudinal sampling scheme of the grids: either 1 for `nlon` = `nlat` or 2 for `nlon` = 2 * `nlat`. |
| `epoch` | The epoch time of the gravity model. |


## Class methods

| Method | Description |
| ------ | ----------- |
| `plot()` | Plot the geoid.|
| `to_xarray()` | Return the gridded data as an xarray DataArray.|
| `to_netcdf()` | Return the gridded data as a netcdf formatted file or object.|
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the SHGravGrid instance. |
