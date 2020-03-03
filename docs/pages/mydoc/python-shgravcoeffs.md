---
title: "SHGravCoeffs class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shgravcoeffs.html
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
| SHGravRealCoeffs | Real gravitational potential spherical harmonic coefficient class. |

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = SHGravCoeffs.from_array()` | Initialize using coefficients from an array. |
| `x = SHGravCoeffs.from_random()` | Initialize using random coefficients with a prescribed power spectrum. |
| `x = SHGravCoeffs.from_zeros()` | Initialize with coefficients set to zero. |
| `x = SHGravCoeffs.from_file()` | Initialize using coefficients from a file. |
| `x = SHGravCoeffs.from_netcdf()` | Initialize using coefficients from a netcdf file. |
| `x = SHGravCoeffs.from_shape()` | Initialize using the gravitational potential predicted from surface relief. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `lmax` | The maximum spherical harmonic degree of the coefficients. |
| `coeffs` | The raw coefficients with the specified normalization and phase conventions. |
| `gm` | The gravitational constant times the mass that is associated with the gravitational potential coefficients. |
| `r0` | The reference radius of the gravitational potential coefficients. |
| `omega` | The angular rotation rate of the body. |
| `normalization` | The normalization of the coefficients: `'4pi'`, `'ortho'`, `'schmidt'`, or `'unnorm'`.|
| `csphase` | Defines whether the Condon-Shortley phase is used (`1`) or not (`-1`). |
| `mask` | A boolean mask that is `True` for the permissible values of degree `l` and order `m`. |
| `kind` | The coefficient data type (only `'real'` is permissible). |
| `header` | A list of values from the header line of the input file used to initialize the class. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `degrees()` | Return an array listing the spherical harmonic degrees from `0` to `lmax`. |
| `spectrum()` | Return the spectrum of the function. |
| `set_omega()` | Set the angular rotation rate of the body. |
| `set_coeffs()` | Set coefficients in-place to specified values.|
| `change_ref()` | Return a new class instance referenced to a different gm, r0, or omega. |
| `rotate()` | Rotate the coordinate system used to express the spherical harmonics coefficients and return a new class instance.|
| `convert()` | Return a new class instance using a different normalization convention. |
| `pad()` | Return a new class instance that is zero padded or truncated to a different `lmax`. |
| `expand()` | Calculate the three vector components of the gravity field, the total field, and the gravitational potential, and return an SHGravGrid class instance. |
| `tensor()` | Calculate the 9 components of the gravity tensor and return an SHGravTensor class instance. |
| `geoid()` | Calculate the height of the geoid and return an SHGeoid class instance. |
| `plot_spectrum()` | Plot the spectrum as a function of spherical harmonic degree. |
| `plot_spectrum2d()` | Plot the spectrum of all spherical-harmonic coefficients. |
| `to_array()` | Return an array of spherical harmonics coefficients with a different normalization convention. |
| `to_file()` | Save raw spherical harmonic coefficients to a text or binary file. |
| `to_netcdf()` | Return the coefficient data as a netcdf formatted file or object. |
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the SHGravCoeffs instance.|