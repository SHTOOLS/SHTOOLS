---
title: "Slepian class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-slepian.html
summary: 
toc: true
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:65%;
}
</style>

## Subclasses

| Subclass name | Description |
| ------------- | ----------- |
| SlepianCap | Class for Slepian functions concentrated within a spherical cap.|
| SlepianMask | Class for Slepian functions concentrated in an arbitrary domain.|

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = Slepian.from_cap()` | Construct Slepian functions concentrated within a spherical cap. |
| `x = Slepian.from_mask()` | Construct Slepian functions concentrated within an arbitrary region. |

## Attributes

| Attribute | Description |
| --------- | ----------- |
| `kind` | Either `'cap'` or `'mask'`.|
| `coeffs` | Array of spherical harmonic coefficients of the rotated spherical-cap Slepian functions. |
| `shannon` | The Shannon number, which approximates the number of well localized functions. |
| `area` | Area of the concentration domain, in radians. |
| `eigenvalues` | Concentration factors of the Slepian functions. |
| `orders` | The angular orders for each of the spherical-cap Slepian functions. |
| `lmax` | Spherical harmonic bandwidth of the Slepian functions. |
| `theta` | Angular radius of the spherical-cap localization domain (default in degrees).|
| `theta_degrees` | `True` (default) if theta is in degrees. |
| `nmax` | Number of Slepian functions. Default is `(lmax+1)**2`. |
| `nrot` | The number of best-concentrated spherical cap Slepian functions that were rotated and whose coefficients are stored in `coeffs`. |
| `clat`, `clon` | Latitude and longitude of the center of the rotated spherical-cap Slepian functions (default in degrees). |
| `coord_degrees` | `True` (default) if `clat` and `clon` are in degrees.|
| `slepian_degrees` | Boolean or int array defining which spherical harmonic degrees were used to construct the Slepian functions. |

## Methods

| Method | Description |
| ------ | ----------- |
| `expand()` | Expand the input function in Slepian functions.|
| `to_array()` | Return an array of the spherical harmonic coefficients for function `alpha`, where `alpha`=`0` is the best concentrated, optionally using a different normalization convention. |
| `to_shcoeffs()` | Return the spherical harmonic coefficients of function `alpha`, where `alpha`=`0` is the best concentrated, as a new SHCoeffs class instance, optionally using a different normalization convention.|
| `to_shgrid()` | Return as a new SHGrid instance a grid of function `alpha`, where `alpha`=`0` is the best concentrated. |
| `number_concentrated()` | Return the number of Slepian functions that have concentration factors greater or equal to a specified value. |
| `degrees()` | Return an array containing the spherical harmonic degrees of the Slepian functions, from `0` to `lmax`. |
| `spectra()` | Return the spectra of one or more Slepian functions.|
| `rotate()` | Rotate the spherical cap Slepian functions, originally located at the North pole, to `clat` and `clon` and save the spherical harmonic coefficients in `coeffs`.|
| `variance()` | Calculate the theoretical variance of the power of a function expanded in spherical-cap Slepian functions. |
| `copy()` | Return a copy of the class instance. |
| `plot()` | Plot the best concentrated Slepian functions using a simple cylindrical projection. |
| `plot_spectra()` | Plot the spectra of the best concentrated Slepian functions. |
| `info()` | Print a summary of the data stored in the Slepian instance. |