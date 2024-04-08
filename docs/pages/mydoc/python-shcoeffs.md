---
title: "SHCoeffs class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shcoeffs.html
summary: Classes for spherical harmonic coefficients.
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:60%;
}
</style>

## Subclasses

| Subclass name | Description |
| ------------- | ----------- |
| SHRealCoeffs | Real spherical harmonic coefficients class. |
| SHComplexCoeffs | Complex spherical harmonic coefficients class. |

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = SHCoeffs.from_array()` | Initialize using coefficients from an array. |
| `x = SHCoeffs.from_random()` | Initialize using random coefficients with a prescribed power spectrum. |
| `x = SHCoeffs.from_zeros()` | Initialize with coefficients set to zero. |
| `x = SHCoeffs.from_file()` | Initialize using coefficients from a file. |
| `x = SHCoeffs.from_netcdf()` | Initialize using coefficients from a netcdf file. |
| `x = SHCoeffs.from_cap()` | Initialize using coefficients of a spherical cap. |
| `x = SHCoeffs.from_least_squares()` | Initialize using coefficients from a least squares inversion. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `lmax` | The maximum spherical harmonic degree of the coefficients. |
| `coeffs` | The raw coefficients with the specified normalization and phase conventions. |
| `errors` | The uncertainties of the spherical harmonic coefficients. |
| `error_kind` | An arbitrary string describing the kind of errors, such as `None`, `'unspecified'`, `'calibrated'` or `'formal'`. |
| `normalization` | The normalization of the coefficients: `'4pi'`, `'ortho'`, `'schmidt'`, or `'unnorm'`. |
| `csphase` | Defines whether the Condon-Shortley phase is used (`1`) or not (`-1`). |
| `mask` | A boolean mask that is `True` for the permissible values of degree `l` and order `m`. |
| `kind` | The coefficient data type: either `'complex'` or `'real'`. |
| `units` | The units of the spherical harmonic coefficients. |
| `header` | A list of values from the header line of the input file used to initialize the class. |
| `header2` | A list of values from the second header line of the input file used to initialize the class. |
| `chi2` | The (weighted) residual sum of squares misfit when initializing from a least squares inversion. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `degrees()` | Return an array listing the spherical harmonic degrees from `0` to `lmax`. |
| `spectrum()` | Return the spectrum of the function. |
| `cross_spectrum()` | Return the cross-spectrum of two functions. |
| `admittance()` | Return the admittance with another function. |
| `correlation()` | Return the spectral correlation with another function. |
| `admitcorr()` | Return the admittance and spectral correlation with another function. |
| `volume()` | Calculate the volume of the body. |
| `centroid()` | Calculate the centroid of the body. |
| `set_coeffs()` | Set coefficients in-place to specified values. |
| `rotate()` | Rotate the coordinate system used to express the spherical harmonics coefficients and return a new class instance. |
| `convert()` | Return a new class instance using a different normalization convention. |
| `pad()` | Return a new class instance that is zero padded or truncated to a different `lmax`. |
| `expand()` | Evaluate the coefficients either on a spherical grid and return an SHGrid class instance, or for a list of latitude and longitude coordinates. |
| `gradient()` | Compute the horizontal gradient of the function and return an SHGradient class instance. |
| `plot_spectrum()` | Plot the spectrum as a function of spherical harmonic degree. |
| `plot_cross_spectrum()` | Plot the cross-spectrum of two functions. |
| `plot_spectrum2d()` | Plot the spectrum of all spherical-harmonic coefficients. |
| `plot_cross_spectrum2d()` | Plot the cross-spectrum of all spherical-harmonic coefficients. |
| `plot_admittance()` | Plot the admittance with another function. |
| `plot_correlation()` | Plot the correlation with another function. |
| `plot_admitcorr()` | Plot the admittance and/or correlation with another function. |
| `to_array()` | Return an array of spherical harmonics coefficients with a different normalization convention. |
| `to_file()` | Save raw spherical harmonic coefficients to a text or binary file. |
| `to_netcdf()` | Return the coefficient data as a netcdf formatted file or object. |
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the SHCoeffs instance. |
