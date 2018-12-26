---
title: "SHGrid class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shgrid.html
summary: 
toc: true
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

## Subclasses

| Subclass name | Description |
| ------------- | ----------- |
| DHRealGrid | Class for real *Driscoll and Healy* (1994) sampled grids.|
| DHComplexGrid | Class for complex *Driscoll and Healy* (1994) sampled grids. |
| GLQRealGrid | Class for real Gauss-Legendre quadrature sampled grids.| 
| GLQComplexGrid | Class for complex Gauss-Legendre quadrature sampled grids.|

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = SHGrid.from_array()` | Initialize using an array. |
| `x = SHGrid.from_file()` | Initialize using an array from a file. |
| `x = SHGrid.from_zeros()` | Initialize using an array of zeros. |
| `x = SHGrid.from_cap()` | Initialize using a rotated spherical cap. |


## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `data` | Array of the gridded data. |
| `nlat`, `nlon` | The number of latitude and longitude bands in the grid.|
| `lmax` | The maximum spherical harmonic degree that can be resolved by the grid sampling. |
| `sampling` | For Driscoll and Healy grids, the longitudinal sampling of the grid. Either 1 for `nlon` = `nlat` or 2 for `nlon` = 2 * `nlat` |
| `kind` | Either `'complex'` or `'real'` for the data type. |
| `grid` | Either `'DH'` or `'GLQ'` for Driscoll and Healy grids or Gauss-Legendre quadrature grids. |
| `zeros` | The $$\cos(\theta)$$ nodes used with Gauss-Legendre quadrature grids. Default is `None`.|
| `weights` | The latitudinal weights used with Gauss-Legendre quadrature grids. Default is `None`. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `to_file()` | Save raw gridded data to a text or binary file. |
| `to_array()` | Return a numpy array of the gridded data. |
| `lats()` | Return a vector containing the latitudes of each row of the gridded data. |
| `lons()` | Return a vector containing the longitudes of each column of the gridded data. |
| `expand()` | Expand the grid into spherical harmonics. |
| `min()` | Return the minimum value of data. |
| `max()` | Return the maximum value of data. |
| `copy()` | Return a copy of the class instance. |
| `plot()` | Plot the raw data using a simple cylindrical projection. |
| `plot3d()` | Plot the raw data on a 3d sphere. |
| `info()` | Print a summary of the data stored in the SHGrid instance. |
