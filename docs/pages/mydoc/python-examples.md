---
title: "Tutorials & guides"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-examples.html
summary: The easiest way to start learning about pyshtools is to see it in action in a jupyter notebook. Start with the tutorials, and then move on to the guides.
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

## Tutorials

The tutorials are designed for python users who are encountering pyshtools for the first time. Each jupyter notebook takes the user by the hand and shows how to perform the most basic operations. Each tutorial should take about 15 minutes to read, and at the end of each you will be able to do simple tasks, like make grids of data from spherical harmonic coefficients, and project these in map form.

| Tutorial | Description |
| ------------- | ----------- |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/grids-and-coefficients.ipynb" target="_blank" rel="noopener">Spherical harmonic coefficients and grids</a> | Learn how to transform spherical harmonic coefficients into maps, maps into spherical harmonic coefficients, and how to plot the power spectrum. |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/localized-spectral-analysis.ipynb" target="_blank" rel="noopener">Localization windows and spectral analysis</a> | Learn how to obtain the power spectrum of a function, localized to any region on the sphere. |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/gravity-and-magnetic-fields.ipynb" target="_blank" rel="noopener">Gravity and magnetic fields</a> | Learn how to read gravity spherical harmonic coefficients from a file, and to make maps of the geoid and free-air gravity. |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/plotting-maps.ipynb" target="_blank" rel="noopener">Plotting maps</a> | Learn how to make publication quality images using geographical projections. |


## Guides

These guides assume that the user already has a basic understanding of how pyshtools works. Each notebook describes in a higher level of detail how to use advanced features of pyshtools that often arise in scientific analyses.

| Guide | Description |
| ------------- | ----------- |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/low-level-spherical-harmonic-analyses.ipynb" target="_blank" rel="noopener">Low-level spherical harmonic analyses</a> | Learn how to do spherical harmonic analyses using low-level functions (without SHCoeffs and SHGrid classes). |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/advanced-shcoeffs-and-shgrid-usage.ipynb" target="_blank" rel="noopener">Advanced usage of SHCoeffs and SHGrid</a> | Learn advanced features of the SHCoeffs and SHGrid class interfaces. |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/spherical-harmonic-normalizations.ipynb" target="_blank" rel="noopener">Spherical harmonic normalizations</a> | Learn more about spherical harmonic normalizations and Parseval's theorem. |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/advanced-localized-spectral-analysis.ipynb" target="_blank" rel="noopener">Advanced localized spectral analysis</a> | Learn more about performing localized spectral analyses on the sphere. |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/advanced-shwindow-usage.ipynb" target="_blank" rel="noopener">Advanced usage of SHWindow</a> | Learn advanced features of the SHWindow class interface. |
| <a href="https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/develop/examples/notebooks/3d-plots.ipynb" target="_blank" rel="noopener">3D plots</a> | Learn how to make 3-dimensional plots of gridded data. |
