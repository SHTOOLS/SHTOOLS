---
title: "Example programs"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: fortran-examples.html
summary: If you want to learn how to incorporate shtools routines in your fortran programs, the following example programs are a good starting point to see shtools in action.
toc: false
folder: fortran
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

{% include note.html content="In order to access the fortran example programs and example datasets, it will be necessary to dowload the entire shtools repo from [GitHub](https://github.com/SHTOOLS/SHTOOLS/), install with *brew* using the option `--with-examples`, or install with *macports*. If the entire repo was downloaded, the example programs will be found in the folder `examples/fortran`. If shtools was instead installed with *brew* or *macports*, they will be found in `/usr/local/share/shtools/examples/fortran/` or `/opt/local/share/shtools/examples/fortran/`, respectively." %}


| Folder | Description |
| ------------- | ----------- |
| `SHCilmPlus/` | Demonstration of how to expand spherical harmonic files into gridded maps using the GLQ routines, and how to compute the gravity field resulting from finite amplitude surface relief. |
| `SHExpandDH/` | Demonstration of how to expand a grid that is equally sampled in latitude and longitude into spherical harmonics using the sampling theorem of *Driscoll and Healy* (1994). |
| `SHExpandLSQ/` | Demonstration of how to expand a set of irregularly sampled data points in latitude and longitude into spherical harmonics by use of a least squares inversion. |
| `SHMag/` | Demonstration of how to expand scalar magnetic potential spherical harmonic coefficients into their three vector components and total field. |
| `MarsCrustalThickness`/ | Demonstration of how to compute a crustal thickness map of Mars. |
| `SHRotate/` | Demonstration of how to determine the spherical harmonic coefficients for a body that is rotated with respect to its initial configuration. |
| `SHLocalizedAdmitCorr/` | Demonstration of how to calculate localized admittance and correlation spectra for a given set of gravity and topography spherical harmonic coefficients. |
| `TimingAccuracy/` | Test programs that calculate the time required to perform the GLQ and DH spherical harmonic transforms and reconstructions and the accuracy of these operations. |
