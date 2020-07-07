---
title: MakeGridDHC (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: makegriddhc.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Create a 2D complex map from a set of complex spherical harmonic coefficients that conforms with Driscoll and Healy's (1994) sampling theorem.

## Usage

call MakeGridDHC (`griddh`, `n`, `cilm`, `lmax`, `norm`, `sampling`, `csphase`, `lmax_calc`, `extend`, `exitstatus`)

## Parameters

`griddh` : output, complex(dp), dimension (nlat, nlong)
:   A 2D complex map of the input spherical harmonic coefficients `cilm` that conforms to the sampling theorem of Driscoll and Healy (1994). If `sampling` is 1, the grid is equally sampled and is dimensioned as (`n` by `n`), where `n` is `2lmax+2`. If sampling is 2, the grid is equally spaced and is dimensioned as (`n` by 2`n`). The first latitudinal band of the grid corresponds to 90 N, the latitudinal sampling interval is 180/`n` degrees, and the default behavior is to exclude the latitudinal band for 90 S. The first longitudinal band of the grid is 0 E, by default the longitudinal band for 360 E is not included, and the longitudinal sampling interval is 360/`n` for an equally sampled and 180/`n` for an equally spaced grid, respectively. If `extend` is 1, the longitudinal band for 360 E and the latitudinal band for 90 S will be included, which increases each of the dimensions of the grid by 1.

`n` : output, integer
:   The number of samples of the grid, which is equal to `2lmax+2`.

`cilm` : input, complex(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The complex spherical harmonic coefficients of the function.  The first index specifies the coefficient corresponding to the positive and negative order of `m`, respectively, with `Clm=cilm(1,l+1,m+1)` and `Cl,-m=cilm(2,l+1,m+1)`.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the function. This determines the number of samples `n`.

`norm` : input, optional, integer, default = 1
:   1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`sampling` : input, optional, integer, default = 1
:   If 1 (default) the input grid is equally sampled (`n` by `n`). If 2, the grid is equally spaced (`n` by 2`n`).

`csphase` : input, optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`lmax_calc` : input, optional, integer, default = `lmax`
:   The maximum spherical harmonic degree used in evaluating the function. This must be less than or equal to `lmax`.

`extend` : input, optional, integer, default = 0
:   If 1, compute the longitudinal band for 360 E and the latitudinal band for 90 S. This increases each of the dimensions of `griddh` by 1.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`MakeGridDHC` will create a 2-dimensional complex map equally sampled (`n` by `n`) or equally spaced (`n` by `2n`) in latitude and longitude from a set of input complex spherical harmonic coefficients, where N is `2lmax+2`. This grid conforms with the sampling theorem of Driscoll and Healy (1994) and this routine is the inverse of `SHExpandDHC`. The function is evaluated at each longitudinal band by inverse Fourier transforming the exponential terms for each degree `l`, and then summing over all degrees. When evaluating the function, the maximum spherical harmonic degree that is considered is the minimum of `lmax`, the size of `cilm-1`, or `lmax_calc` (if specified).

The default is to use an input grid that is equally sampled (`n` by `n`), but this can be changed to use an equally spaced grid (`n` by 2`n`) by the optional argument `sampling`. The redundant longitudinal band for 360 E and the latitudinal band for 90 S are excluded by default, but these can be computed by specifying the optional argument `extend`. The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

The normalized legendre functions are calculated using the scaling algorithm of Holmes and Featherstone (2002), which are accurate to about degree 2800. The unnormalized functions are accurate only to about degree 15.

## References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, J. Geodesy, 76, 279-299, 2002.

## See also

[shexpanddhc](shexpanddhc.html), [makegriddh](makegriddh.html), [shexpanddh](shexpanddh.html), [makegridglq](makegridglq.html), [shexpandglq](shexpandglq.html), [makegridglqc](makegridglqc.html), [shexpandglqc](shexpandglqc.html), [shexpandlsq](shexpandlsq.html)
