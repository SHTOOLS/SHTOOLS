---
title: "Spherical harmonic transforms"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: spherical-harmonic-transforms.html
summary: 
toc: true
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

| Routine name | Description |
| ------------ | ----------- |
| [SHExpandDH](shexpanddh.html) | Expand an equally sampled or equally spaced map into spherical harmonics using *Driscoll and Healy*'s (1994) sampling theorem. |
| [MakeGridDH](makegriddh.html) | Create a 2D map from a set of spherical harmonic coefficients that conforms with *Driscoll and Healy*'s (1994) sampling theorem. |
| [SHExpandDHC](shexpanddhc.html) | Expand an equally sampled or equally spaced complex map into complex spherical harmonics using *Driscoll and Healy*'s (1994) sampling theorem. |
| [MakeGridDHC](makegriddhc.html) | Create a 2D complex map from a set of complex spherical harmonic coefficients that conforms with *Driscoll and Healy*'s (1994) sampling theorem. |

## Gauss-Legendre quadrature grids

| Routine name | Description |
| ------------ | ----------- |
| [SHGLQ](shglq.html) | Precompute weights, nodes, and associated Legendre functions used in the GLQ-based spherical harmonics routines. |
| [SHExpandGLQ](shexpandglq.html) | Expand a 2D map sampled on the Gauss-Legendre quadrature nodes into spherical harmonics. |
| [MakeGridGLQ](makegridglq.html) | Create a 2D map from a set of spherical harmonic coefficients sampled on a the Gauss-Legendre quadrature nodes. |
| [SHExpandGLQC](shexpandglqc.html) | Expand a 2D complex map sampled on the Gauss-Legendre quadrature nodes into complex spherical harmonics. |
| [MakeGridGLQC](makegridglqc.html) | Create a 2D complex map from a set of complex spherical harmonic coefficients sampled on a the Gauss-Legendre quadrature nodes. |
| [GLQGridCoord](glqgridcoord.html) | Compute the latitude and longitude coordinates used in Gauss-Legendre quadrature grids. |

## Other routines

| Routine name | Description |
| ------------ | ----------- |
| [SHExpandLSQ](shexpandlsq.html) | Expand a set of irregularly sampled data points into spherical harmonics using a (weighted) least squares inversion. |
| [MakeGrid2D](makegrid2d.html) | Create a 2D cylindrical map with arbitrary grid spacing from a set of spherical harmonic coefficients. |
| [MakeGridPoint](makegridpoint.html) | Evaluate a real function expressed in real spherical harmonics at a single point. |
| [MakeGridPointC](makegridpointc.html) | Evaluate a complex function expressed in complex spherical harmonics at a single point. |
| [SHMultiply](shmultiply.html) | Multiply two functions and determine the spherical harmonic coefficients of the resulting function. |

