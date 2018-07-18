---
title: "SHCoeffs class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shcoeffs.html
summary: 
toc: true
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
| SHRealCoeffs | Real spherical harmonic coefficient class. |
| SHComplexCoeffs | Complex spherical harmonic coefficient class. |

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = SHCoeffs.from_array()` | Initialize using coefficients from an array. |
| `x = SHCoeffs.from_random()` | Initialize using random coefficients with a prescribed power spectrum. |
| `x = SHCoeffs.from_zeros()` | Initialize with coefficients set to zero. |
| `x = SHCoeffs.from_file()` | Initialize using coefficients from a file. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `lmax` | The maximum spherical harmonic degree of the coefficients. |
| `coeffs` | The raw coefficients with the specified normalization and phase conventions. |
| `normalization` | The normalization of the coefficients: `'4pi'`, `'ortho'`, or `'schmidt'`.|
| `csphase` | Defines whether the Condon-Shortley phase is used (`1`) or not (`-1`). |
| `mask` | A boolean mask that is `True` for the permissible values of degree `l` and order `m`. |
| `kind` | The coefficient data type: either `'complex'` or `'real'`. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `to_file()` | Save raw spherical harmonic coefficients to a text or binary file. |
| `to_array()` | Return an array of spherical harmonics coefficients with a different normalization convention. |
| `degrees()` | Return an array listing the spherical harmonic degrees from `0` to `lmax`. |
| `spectrum()` | Return the spectrum of the function.|
| `set_coeffs()` | Set coefficients in-place to specified values.|
| `rotate()` | Rotate the coordinate system used to express the spherical harmonics coefficients and return a new class instance.|
| `convert()` | Return a new class instance using a different normalization convention. |
| `pad()` | Return a new class instance that is zero padded or truncated to a different `lmax`.|
| `expand()` | Evaluate the coefficients either on a spherical grid and return an SHGrid class instance, or for a list of latitude and longitude coordinates.| 
| `copy()` | Return a copy of the class instance. |
| `plot_spectrum()` | Plot the spectrum as a function of spherical harmonic degree. |
| `plot_spectrum2d()` | Plot the spectrum of all spherical-harmonic coefficients. |
| `info()` | Print a summary of the data stored in the SHCoeffs instance.|
