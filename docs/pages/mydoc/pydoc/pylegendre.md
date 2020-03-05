---
title: legendre (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pylegendre.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute all the associated Legendre functions up to a maximum degree and order.

## Usage

`plm` = legendre (`lmax`, `z`, [`normalization`, `csphase`, `cnorm`, `packed`])

## Returns

`plm` : float, dimension (`lmax`+1, `lmax`+1) or ((`lmax`+1)\*(`lmax`+2)/2)
:   An array of associated Legendre functions, plm[l, m], where `l` and `m` are the degree and order, respectively. If `packed` is True, the array is 1-dimensional with the index corresponding to `l\*(l+1)/2+m`.

## Parameters

`lmax` : integer
:   The maximum degree of the associated Legendre functions to be computed.

`z` : float
:   The argument of the associated Legendre functions.

`normalization` : str, optional, default = '4pi'
:   '4pi', 'ortho', 'schmidt', or 'unnorm' for use with geodesy 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized spherical harmonic functions, respectively.

`csphase` : optional, integer, default = 1
:   If 1 (default), the Condon-Shortley phase will be excluded. If -1, the Condon-Shortley phase of (-1)^m will be appended to the associated Legendre functions.

`cnorm` : optional, integer, default = 0
:   If 1, the complex normalization of the associated Legendre functions will be used. The default is to use the real normalization.

`packed` : optional, bool, default = False
:   If True, return a 1-dimensional packed array with the index corresponding to `l\*(l+1)/2+m`, where `l` and `m` are respectively the degree and order.

## Description

`legendre` will calculate all of the associated Legendre functions up to degree `lmax` for a given argument. The Legendre functions are used typically as a part of the spherical harmonic functions, and three parameters determine how they are defined. `normalization` can be either '4pi' (default), 'ortho', 'schmidt', or 'unnorm' for use with 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized spherical harmonic functions, respectively. `csphase` determines whether to include or exclude (default) the Condon-Shortley phase factor. `cnorm` determines whether to normalize the Legendre functions for use with real (default) or complex spherical harmonic functions.

By default, the routine will return a 2-dimensional array, p[l, m]. If the optional parameter `packed` is set to True, the output will instead be a 1-dimensional array where the indices correspond to `l\*(l+1)/2+m`. The Legendre functions are calculated using the standard three-term recursion formula, and in order to prevent overflows, the scaling approach of Holmes and Featherstone (2002) is utilized. The resulting functions are accurate to about degree 2800. See Wieczorek and Meschede (2018) for exact definitions on how the Legendre functions are defined.

## References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, J. Geodesy, 76, 279-299, doi:10.1007/s00190-002-0216-2, 2002.

Wieczorek, M. A., and M. Meschede. SHTools - Tools for working with spherical harmonics, Geochem., Geophys., Geosyst., 19, 2574-2592, doi:10.1029/2018GC007529, 2018.

## See also

[pylegendre_lm](pylegendre_lm.html), [pyplmbar](pyplmbar.html), [pyplmon](pyplmon.html), [pyplmschmidt](pyplmschmidt.html), [pyplegendrea](pyplegendrea.html), [pyplmindex](pyplmindex.html)

