---
title: "SHWindow class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shwindow.html
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
| SHWindowCap | Class for windows concentrated within a spherical cap. |
| SHWindowMask | Class for windows concentrated in an arbitrary domain. |

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = SHWindow.from_cap()` | Construct windows concentrated within a spherical cap. |
| `x = SHWindow.from_mask()` | Construct windows concentrated within an arbitrary region. |

## Attributes

| Attribute | Description |
| --------- | ----------- |
| `kind` | Either `'cap'` or `'mask'`. |
| `coeffs` | Array of spherical harmonic coefficients of the rotated spherical-cap localization windows. |
| `shannon` | The Shannon number, which approximates the number of well localized windows. |
| `area` | Area of the concentration domain, in radians. |
| `eigenvalues` | Concentration factors of the localization windows. |
| `orders` | The angular orders for each of the spherical-cap localization windows. |
| `weights` | Taper weights used with the multitaper spectral analyses. Default is `None`. |
| `lmax` | Spherical harmonic bandwidth of the localization windows. |
| `theta` | Angular radius of the spherical-cap localization domain (default in degrees). |
| `theta_degrees` | `True` (default) if theta is in degrees. |
| `nwin` | Number of localization windows. Default is `(lmax+1)**2`. |
| `nwinrot` | The number of best-concentrated spherical cap windows that were rotated and whose coefficients are stored in `coeffs`. |
| `clat`, `clon` | Latitude and longitude of the center of the rotated spherical-cap localization windows (default in degrees). |
| `coord_degrees` | `True` (default) if `clat` and `clon` are in degrees. |

## Methods

| Method | Description |
| ------ | ----------- |
| `to_array()` | Return an array of the spherical harmonic coefficients for taper `i`, where `i`=`0` is the best concentrated, optionally using a different normalization convention. |
| `to_shcoeffs()` | Return the spherical harmonic coefficients of taper `i`, where `i`=`0` is the best concentrated, as a new SHCoeffs class instance, optionally using a different normalization convention. |
| `to_shgrid()` | Return as a new SHGrid instance a grid of taper `i`, where `i`=`0`  is the best concentrated window. |
| `number_concentrated()` | Return the number of localization windows that have concentration factors greater or equal to a specified value. |
| `degrees()` | Return an array containing the spherical harmonic degrees of the localization windows, from `0` to `lmax`. |
| `spectra()` | Return the spectra of one or more localization windows. |
| `rotate()` | Rotate the spherical cap tapers, originally located at the north pole, to `clat` and `clon` and save the spherical harmonic coefficients in `coeffs`. |
| `coupling_matrix()` | Return the coupling matrix of the first `nwin`. |
| `biased_spectrum()` | Calculate the multitaper (cross-) spectrum expectation of a localized function. |
| `multitaper_spectrum()` | Return the multitaper spectrum estimate and uncertainty for the input SHCoeffs class instance. |
| `multitaper_cross_spectrum()` | Return the multitaper cross-spectrum estimate and uncertainty for two input SHCoeffs class instances. |
| `copy()` | Return a copy of the class instance. |
| `plot_windows()` | Plot the best concentrated localization windows using a simple cylindrical projection. |
| `plot_spectra()` | Plot the spectra of the best concentrated localization windows. |
| `plot_coupling_matrix()` | Plot the multitaper coupling matrix.|
| `info()` | Print a summary of the data stored in the SHWindow instance. |