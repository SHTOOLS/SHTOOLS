---
title: "Python examples"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-examples.html
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

## Notebooks

| Notebook name | Description |
| ------------- | ----------- |
| <a href="pages/mydoc/notebooks/Introduction-1.html" target="_blank" rel="noopener">Introduction 1</a> | Grids and and Spherical Harmonic Coefficients. |
| <a href="pages/mydoc/notebooks/Introduction-2.html" target="_blank" rel="noopener">Introduction 2</a> | Localization windows and spectral analysis. |
| <a href="pages/mydoc/notebooks/tutorial_1.html" target="_blank" rel="noopener">Tutorial 1</a> | Simple spherical harmonic analyses. |
| <a href="pages/mydoc/notebooks/tutorial_2.html" target="_blank" rel="noopener">Tutorial 2</a> | Localized spectral analysis on the sphere. |
| <a href="pages/mydoc/notebooks/tutorial_3.html" target="_blank" rel="noopener">Tutorial 3</a> | The SHTOOLS class interface. |
| <a href="pages/mydoc/notebooks/tutorial_4.html" target="_blank" rel="noopener">Tutorial 4</a> | Spherical harmonic normalizations and Parseval's theorem. |
| <a href="pages/mydoc/notebooks/tutorial_5.html" target="_blank" rel="noopener">Tutorial 5</a> | Multitaper spectral analysis class interface. |
| <a href="pages/mydoc/notebooks/tutorial_6.html" target="_blank" rel="noopener">Tutorial 6</a> | 3D plots of gridded data.|

## Test programs

A variety of test programs can be found in the folders in `examples/python`.

| Folder directory | Description |
| ------------- | ----------- |
| `ClassInterface/` | Test the python class interfaces. |
| `TestLegendre/` | Test and plot the Legendre functions. |
| `IOStorageConversions/` | Read coefficients from a file and test conversions between real and complex coefficients. |
| `GlobalSpectralAnalysis/` | Test functions to compute different power spectra from real and complex coefficients. |
| `LocalizedSpectralAnalysis/` | Test the coupling matrix, localized spectral analysis, and bias routines. |
| `GravMag/` | Test the gravity and magnetics routines, and compute the crustal thickness of Mars.|
| `TimingAccuracy/` | Perform timing and accuracy tests using real and complex coefficients, with *Driscoll and Healy* (1994) and Gauss-Lengendre quadrature grids.|
| `Other/` | Test a variety of other routines.|
