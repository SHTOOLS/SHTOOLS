---
title: "SHMagGrid class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shmaggrid.html
summary: A class for global gridded magnetic field data.
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
| `x = SHMagCoeffs.expand()` | Initialize using an SHMagCoeffs class instance. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `rad` | SHGrid class instance of the radial component of the magnetic field evaluated on an ellipsoid. |
| `theta` | SHGrid class instance of the theta component of the magnetic field evaluated on an ellipsoid. |
| `phi` | SHGrid class instance of the phi component of the magnetic field evaluated on an ellipsoid. |
| `total` | SHGrid class instance of the total magnetic field evaluated on an ellipsoid. |
| `pot` | SHGrid class instance of the magnetic potential evaluated on an ellipsoid. |
| `a` | Semimajor axis of the reference ellipsoid. |
| `f` | Flattening of the reference ellipsoid, f = (a - b) / a. |
| `lmax` | The maximum spherical harmonic degree resolvable by the grids. |
| `lmax_calc` | The maximum spherical harmonic degree of the magnetic potential used in creating the grids. |
| `units` | The units of the gridded magnetic field data. |
| `pot_units` | The units of the gridded magnetic potential data. |
| `nlat`, `nlon` | The number of latitude and longitude bands in the grids. |
| `sampling` | The longitudinal sampling scheme of the grids: either 1 for `nlon` = `nlat` or 2 for `nlon` = 2 * `nlat`. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `plot()` | Plot all three components of the magnetic field and total magnetic intensity.|
| `plot_rad()` | Plot the radial component of the magnetic field. |
| `plot_theta()` | Plot the theta component of the magnetic field. |
| `plot_phi()` | Plot the phi component of the magnetic field. |
| `plot_total()` | Plot the total magnetic field intensity. |
| `plot_pot()` | Plot the magnetic potential. |
| `to_xarray()` | Return the magnetic field gridded data as an xarray DataSet. |
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the SHMagGrid instance. |
