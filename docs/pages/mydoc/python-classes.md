---
title: "Classes"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-classes.html
summary: 
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>

## Generic classes

| Class name | Description |
| ---------- | ----------- |
| [SHCoeffs](python-shcoeffs.html) | Class for spherical harmonic coefficients.|
| [SHGrid](python-shgrid.html) | Class for global gridded data.|
| [SHWindow](python-shwindow.html) | Class for localized spectral analyses.|
| [Slepian](python-slepian.html) | Class for Slepian functions.|
| [SlepianCoeffs](python-slepiancoeffs.html) | Class for Slepian expansion coefficients.|

## Classes for gravity fields

| [SHGravCoeffs](python-shgravcoeffs.html) | Class for gravitational potential spherical harmonic coefficients.|
| [SHGravGrid](python-shgravgrid.html) | Class for global gridded gravitational field data.|
| [SHGravTensor](python-shtensor.html) | Class for the gravity field tensor and eigenvalues. |
| [SHGeoid](python-shgeoid.html) | Class for the geoid.|

## Classes for magnetic fields

| [SHMagCoeffs](python-shmagcoeffs.html) | Class for magnetic potential spherical harmonic coefficients.|
| [SHMagGrid](python-shmaggrid.html) | Class for global gridded magnetic field data.|
| [SHMagTensor](python-shtensor.html) | Class for the magnetic field tensor and eigenvalues.|
