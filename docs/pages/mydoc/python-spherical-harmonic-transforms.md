---
title: "Spherical harmonic transforms"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-spherical-harmonic-transforms.html
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

## Equally sampled (N&#215;N) and equally spaced (N&#215;2N) grids

| Function name | Description |
| ------------- | ----------- |
| [SHExpandDH](pyshexpanddh.html) | Expand an equally sampled or equally spaced map into spherical harmonics using *Driscoll and Healy*'s (1994) sampling theorem. |
| [MakeGridDH](pymakegriddh.html) | Create a 2D map from a set of spherical harmonic coefficients that conforms with *Driscoll and Healy*'s (1994) sampling theorem. |
| [SHExpandDHC](pyshexpanddhc.html) | Expand an equally sampled or equally spaced complex map into complex spherical harmonics using *Driscoll and Healy*'s (1994) sampling theorem. |
| [MakeGridDHC](pymakegriddhc.html) | Create a 2D complex map from a set of complex spherical harmonic coefficients that conforms with *Driscoll and Healy*'s (1994) sampling theorem. |

## Gauss-Legendre quadrature grids

| Function name | Description |
| ------------- | ----------- |
| [SHGLQ](pyshglq.html) | Precompute the weights and nodes used in the GLQ-based spherical harmonics routines. |
| [SHExpandGLQ](pyshexpandglq.html) | Expand a 2D map sampled on the Gauss-Legendre quadrature nodes into spherical harmonics. |
| [MakeGridGLQ](pymakegridglq.html) | Create a 2D map from a set of spherical harmonic coefficients sampled on a the Gauss-Legendre quadrature nodes. |
| [SHExpandGLQC](pyshexpandglqc.html) | Expand a 2D complex map sampled on the Gauss-Legendre quadrature nodes into complex spherical harmonics. |
| [MakeGridGLQC](pymakegridglqc.html) | Create a 2D complex map from a set of complex spherical harmonic coefficients sampled on a the Gauss-Legendre quadrature nodes. |
| [GLQGridCoord](pyglqgridcoord.html) | Compute the latitude and longitude coordinates used in Gauss-Legendre quadrature grids. |

## Other routines

| Function name | Description |
| ------------- | ----------- |
| [SHExpandLSQ](pyshexpandlsq.html) | Expand a set of irregularly sampled data points into spherical harmonics using a least squares inversion. |
| [SHExpandWLSQ](pyshexpandwlsq.html) | Expand a set of irregularly sampled data points into spherical harmonics using a weighted least squares inversion. |
| [MakeGrid2D](pymakegrid2d.html) | Create a 2D cylindrical map with arbitrary grid spacing from a set of spherical harmonic coefficients. |
| [MakeGridPoint](pymakegridpoint.html) | Evaluate a real function expressed in real spherical harmonics at a single point. |
| [MakeGridPointC](pymakegridpointc.html) | Evaluate a complex function expressed in complex spherical harmonics at a single point. |
| [SHMultiply](pyshmultiply.html) | Multiply two functions and determine the spherical harmonic coefficients of the resulting function. |
| [spharm](pyspharm.html) | Compute all the spherical harmonic functions up to a maximum degree and order. |
| [spharm_lm](pyspharm_lm.html) | Compute the spherical harmonic function for a specific degree l and order m. |
