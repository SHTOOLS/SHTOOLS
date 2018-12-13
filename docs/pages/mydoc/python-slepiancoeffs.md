---
title: "SlepianCoeffs class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-slepiancoeffs.html
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

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = Slepian.expand()` | Determine Slepian expansion coefficients for an input function. |

## Attributes

| Attribute | Description |
| --------- | ----------- |
| `falpha` | Array of the Slepian expansion coefficients. |
| `galpha` | A Slepian class instance that contains the associated Slepian functions. |
| `nmax` | The number of Slepian expansion coefficients. |


## Methods

| Method | Description |
| ------ | ----------- |
| `expand()` | Expand the function on a grid an return an SHGrid class instance. 
| `to_shcoeffs()` | Return the spherical harmonic coefficients of the function as an SHCoeffs class instance.|
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the Slepian instance. |