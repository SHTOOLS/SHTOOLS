---
title: "Fortran examples"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: fortran-examples.html
summary: 
toc: false
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

A variety of test programs can be found in the folders in `examples/fortran`.

| Folder directory | Description |
| ------------- | ----------- |
| `SHCilmPlus/` | Demonstration of how to expand spherical harmonic files into gridded maps using the GLQ routines, and how to compute the gravity field resulting from finite amplitude surface relief. |
| `SHExpandDH/` | Demonstration of how to expand a grid that is equally sampled in latitude and longitude into spherical harmonics using the sampling theorem of *Driscoll and Healy* (1994). |
| `SHExpandLSQ/` | Demonstration of how to expand a set of irregularly sampled data points in latitude and longitude into spherical harmonics by use of a least squares inversion. |
| `SHMag/` | Demonstration of how to expand scalar magnetic potential spherical harmonic coefficients into their three vector components and total field. |
| `MarsCrustalThickness`/ | Demonstration of how to compute a crustal thickness map of Mars. |
| `SHRotate/` | Demonstration of how to determine the spherical harmonic coefficients for a body that is rotated with respect to its initial configuration. |
| `SHLocalizedAdmitCorr/` | Demonstration of how to calculate localized admittance and correlation spectra for a given set of gravity and topography spherical harmonic coefficients. |
| `TimingAccuracy/` | Test programs that calculate the time required to perform the GLQ and DH spherical harmonic transforms and reconstructions and the accuracy of these operations. |
