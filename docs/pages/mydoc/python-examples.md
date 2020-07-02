---
title: "Tutorials & guides"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-examples.html
summary: The easiest way to start learning about pyshtools is to see it in action in a jupyter notebook. Start with the tutorials, and then move on to the guides and example code.
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

{% include note.html content="In order to access the notebooks and example programs, it will be necessary to download or clone the entire shtools repo from [GitHub](https://github.com/SHTOOLS/SHTOOLS/). These files are not included when installing via pip or conda." %}

## Tutorials

The tutorials are designed for python users who are encountering pyshtools for the first time. Each jupyter notebook takes the user by the hand and shows how to perform the most basic operations. Each tutorial should take about 15 minutes to read, and at the end of each you will be able to do simple tasks, like make grids of data from spherical harmonic coefficients, and project these in map form.

| Tutorial | Description |
| ------------- | ----------- |
| <a href="pages/mydoc/notebooks/grids-and-coefficients.html" target="_blank" rel="noopener">Spherical harmonic coefficients and grids</a> | Learn how to transform spherical harmonic coefficients into maps, maps into spherical harmonic coefficients, and how to plot the power spectrum. |
| <a href="pages/mydoc/notebooks/localized-spectral-analysis.html" target="_blank" rel="noopener">Localization windows and spectral analysis</a> | Learn how to obtain the power spectrum of a function, localized to any region on the sphere. |
| <a href="pages/mydoc/notebooks/gravity-and-magnetic-fields.html" target="_blank" rel="noopener">Gravity and magnetic fields</a> | Learn how to read gravity spherical harmonic coefficients from a file, and to make maps of the geoid and free-air gravity. |
| <a href="pages/mydoc/notebooks/plotting-maps.html" target="_blank" rel="noopener">Plotting maps</a> | Learn how to make publication quality images using geographical projections. |


## Guides

These guides assume that the user already has a basic understanding of how pyshtools works. Each notebook describes in a higher level of detail how to use advanced features of pyshtools that often arise in scientific analyses.

| Guide | Description |
| ------------- | ----------- |
| <a href="pages/mydoc/notebooks/low-level-spherical-harmonic-analyses.html" target="_blank" rel="noopener">Low-level spherical harmonic analyses</a> | Learn how to do spherical harmonic analyses using low-level functions (without SHCoeffs and SHGrid classes). |
| <a href="pages/mydoc/notebooks/advanced-shcoeffs-and-shgrid-usage.html" target="_blank" rel="noopener">Advanced usage of SHCoeffs and SHGrid</a> | Learn advanced features of the SHCoeffs and SHGrid class interfaces. |
| <a href="pages/mydoc/notebooks/spherical-harmonic-normalizations.html" target="_blank" rel="noopener">Spherical harmonic normalizations</a> | Learn more about spherical harmonic normalizations and Parseval's theorem. |
| <a href="pages/mydoc/notebooks/advanced-localized-spectral-analysis.html" target="_blank" rel="noopener">Advanced localized spectral analysis</a> | Learn more about performing localized spectral analyses on the sphere. |
| <a href="pages/mydoc/notebooks/advanced-shwindow-usage.html" target="_blank" rel="noopener">Advanced usage of SHWindow</a> | Learn advanced features of the SHWindow class interface. |
| <a href="pages/mydoc/notebooks/3d-plots.html" target="_blank" rel="noopener">3D plots</a> | Learn how to make 3-dimensional plots of gridded data. |


## Test programs

The test programs in `examples/python` demonstrate how pyshtools can be used in real-life situations.

| Folder directory | Description |
| ------------- | ----------- |
| `ClassInterface/` | Test the python class interfaces. |
| `TestLegendre/` | Test and plot the Legendre functions. |
| `IOStorageConversions/` | Read coefficients from a file and test conversions between real and complex coefficients. |
| `GlobalSpectralAnalysis/` | Test functions to compute different power spectra from real and complex coefficients. |
| `LocalizedSpectralAnalysis/` | Test the coupling matrix, localized spectral analysis, and bias routines. |
| `GravMag/` | Test the gravity and magnetics routines, and compute the crustal thickness of Mars.|
| `TimingAccuracy/` | Perform timing and accuracy tests using real and complex coefficients, with *Driscoll and Healy* (1994) and Gauss-Legendre quadrature grids.|
| `Other/` | Test a variety of other routines.|
