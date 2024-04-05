---
title: pysh.expand
keywords: spherical harmonics, spherical harmonic transforms, spherical harmonic expansions, python, pyshtools
sidebar: mydoc_sidebar
permalink: python-spherical-harmonic-transforms.html
summary: This module provides routines for performing spherical harmonic expansions and the construction of grids from spherical harmonic coefficients.
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
| [MakeGradientDH](pymakegradientdh.html) | Compute the gradient of a scalar function and return grids of the two horizontal components that conform with *Driscoll and Healy*'s (1994) sampling theorem. |


## Gauss-Legendre quadrature grids

| Function name | Description |
| ------------- | ----------- |
| [SHGLQ](pyshglq.html) | Precompute the weights and nodes used in the GLQ-based spherical harmonics routines. |
| [SHExpandGLQ](pyshexpandglq.html) | Expand a 2D map sampled on the Gauss-Legendre quadrature nodes into spherical harmonics. |
| [MakeGridGLQ](pymakegridglq.html) | Create a 2D map from a set of spherical harmonic coefficients sampled on a the Gauss-Legendre quadrature nodes. |
| [SHExpandGLQC](pyshexpandglqc.html) | Expand a 2D complex map sampled on the Gauss-Legendre quadrature nodes into complex spherical harmonics. |
| [MakeGridGLQC](pymakegridglqc.html) | Create a 2D complex map from a set of complex spherical harmonic coefficients sampled on a the Gauss-Legendre quadrature nodes. |
| [GLQGridCoord](pyglqgridcoord.html) | Compute the latitude and longitude coordinates used in Gauss-Legendre quadrature grids. |

## Least squares inversion

| Function name | Description |
| ------------- | ----------- |
| [SHExpandLSQ](pyshexpandlsq.html) | Determine the spherical harmonic coefficients of an irregularly sampled function
using a least squares inversion. |
| [SHExpandLSQ_G](pyshexpandlsq_g.html) | Determine the spherical harmonic coefficients of an irregularly sampled function
using a least squares inversion with a precomputed data kernel matrix. |
| [SHExpandWLSQ](pyshexpandwlsq.html) | Determine the spherical harmonic coefficients of an irregularly sampled function
using a weighted least squares inversion. |
| [SHExpandWLSQ_G](pyshexpandwlsq_g.html) | Determine the spherical harmonic coefficients of an irregularly sampled function
using a weighted least squares inversion with a precomputed data kernel matrix. |
| [LSQ_G](pylsq_g.html) | Compute the data kernel matrix G that is used when computing spherical harmonic coefficients by least squares inversion. |

## Other routines

| Function name | Description |
| ------------- | ----------- |
| [MakeGrid2D](pymakegrid2d.html) | Create a 2D cylindrical map with arbitrary grid spacing from a set of spherical harmonic coefficients. |
| [MakeGridPoint](pymakegridpoint.html) | Evaluate a real function expressed in real spherical harmonics at a set of points. |
| [MakeGridPointC](pymakegridpointc.html) | Evaluate a complex function expressed in complex spherical harmonics at a set of points. |
| [SHMultiply](pyshmultiply.html) | Multiply two functions and determine the spherical harmonic coefficients of the resulting function. |
| [spharm](pyspharm.html) | Compute all the spherical harmonic functions up to a maximum degree and order. |
| [spharm_lm](pyspharm_lm.html) | Compute the spherical harmonic function for specific degrees l and orders m. |
