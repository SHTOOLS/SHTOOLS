---
title: legendre_lm()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pylegendre_lm.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute the associated Legendre function for specific degrees and orders.

## Usage

plm = legendre_lm (l, m, z, [normalization, csphase, cnorm])

## Returns

plm : float, ndarray
:   The associated Legendre functions for degree l and order m.

## Parameters

l : integer, array_like
:   The spherical harmonic degree.

m : integer, array_like
:   The spherical harmonic order.

z : float, array_like
:   The argument of the associated Legendre functions.

normalization : str, array_like, optional, default = '4pi'
:   '4pi', 'ortho', 'schmidt', or 'unnorm' for use with geodesy 4pi
    normalized, orthonormalized, Schmidt semi-normalized, or unnormalized
    spherical harmonic functions, respectively.

csphase : integer, array_like, optional, default = 1
:   If 1 (default), the Condon-Shortley phase will be excluded. If -1, the
    Condon-Shortley phase of (-1)^m will be appended to the associated
    Legendre functions.

cnorm : integer, array_like, optional, default = 0
:   If 1, the complex normalization of the associated Legendre functions
    will be used. The default is to use the real normalization.

## Notes

legendre_lm will calculate the associated Legendre function for specific
degrees l and orders m. The Legendre functions are used typically as a part
of the spherical harmonic functions, and three parameters determine how
they are defined. normalization can be either '4pi' (default), 'ortho',
'schmidt', or 'unnorm' for use with 4pi normalized, orthonormalized,
Schmidt semi-normalized, or unnormalized spherical harmonic functions,
respectively. csphase determines whether to include or exclude (default)
the Condon-Shortley phase factor. cnorm determines whether to normalize
the Legendre functions for use with real (default) or complex spherical
harmonic functions.

The Legendre functions are calculated using the standard three-term
recursion formula, and in order to prevent overflows, the scaling approach
of Holmes and Featherstone (2002) is utilized. The resulting functions are
accurate to about degree 2800. See Wieczorek and Meschede (2018) for exact
definitions on how the Legendre functions are defined.

## References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and order
normalised associated Legendre functions, J. Geodesy, 76, 279-299,
doi:10.1007/s00190-002-0216-2, 2002.

Wieczorek, M. A., and M. Meschede. SHTools — Tools for working with
spherical harmonics, Geochem., Geophys., Geosyst., 19, 2574-2592,
doi:10.1029/2018GC007529, 2018.

