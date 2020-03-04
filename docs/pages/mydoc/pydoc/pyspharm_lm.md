---
title: spharm_lm_ (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyspharm_lm.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute the spherical harmonic function for a specific degree and order.

## Usage

`ylm` = spharm (`l`, `m`, `theta`, `phi`, [`normalization`, `kind`, `csphase`, `degrees`])

## Returns

`ylm` : float or complex
:   The spherical harmonic function ylm, where `l` and `m` are the spherical harmonic degree and order, respectively.

## Parameters

`l` : integer
:   The spherical harmonic degree.

`m` : integer
:   The spherical harmonic order.

`theta` : float
:   The colatitude in degrees. Use radians if 'degrees' is set to False.

`phi` : float
:   The longitude in degrees. Use radians if 'degrees' is set to False.

`normalization` : str, optional, default = '4pi'
:   '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized spherical harmonic functions, respectively.

`kind` : str, optional, default = 'real'
:   'real' or 'complex' spherical harmonic coefficients.

`csphase` : optional, integer, default = 1
:   If 1 (default), the Condon-Shortley phase will be excluded. If -1, the Condon-Shortley phase of (-1)^m will be appended to the spherical harmonic functions.

`degrees` : optional, bool, default = True
:   If True, `colat` and `phi` are expressed in degrees.

## Description

`spharm_lm` will calculate the spherical harmonic function for a specific degree `l` and order `m`, and for a given colatitude `theta` and longitude `phi`. Three parameters determine how the spherical harmonic functions are defined. `normalization` can be either '4pi' (default), 'ortho', 'schmidt', or 'unnorm' for 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized spherical harmonic functions, respectively. `kind` can be either 'real' or 'complex', and `csphase` determines whether to include or exclude (default) the Condon-Shortley phase factor.

The spherical harmonic functions are calculated using the standard three-term recursion formula, and in order to prevent overflows, the scaling approach of Holmes and Featherstone (2002) is utilized. The resulting functions are accurate to about degree 2800. See Wieczorek and Meschede (2018) for exact definitions on how the spherical harmonic functions are defined.

## References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, J. Geodesy, 76, 279-299, doi:10.1007/s00190-002-0216-2, 2002.

Wieczorek, M. A., and M. Meschede. SHTools — Tools for working with spherical harmonics, Geochem., Geophys., Geosyst., 19, 2574-2592, doi:10.1029/2018GC007529, 2018.

## See also

[pyspharm](pyspharm.html), [pylegendre_lm](pylegendre_lm.html), [pylegendre](pylegendre.html)
