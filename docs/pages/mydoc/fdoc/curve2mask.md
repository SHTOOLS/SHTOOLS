---
title: Curve2Mask (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: curve2mask.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Given a set of latitude and longitude coordinates representing a closed curve, output a gridded binary mask.

## Usage

call Curve2Mask (`mask_dh`, `n`, `sampling`, `profile`, `nprofile`, `np`, `centralmeridian`, `exitstatus`)

## Parameters

`dh_mask` : output, integer, dimension (`n`, `n`\*`sampling`)
:   A Driscoll and Healy (1994) sampled grid describing the concentration region R. All elements on output will either be 1 (for inside the concentration region) or 0 (for outside R).

`n` : input, integer
:   The number of latitudinal samples in `dh_mask`. The effective spherical harmonic bandwidth of this grid is `L=n/2-1`.

`sampling` : input, integer
:   For 1, `dh_mask` has `n` x `n` samples. For 2, `dh_mask` has `n` x `2n` samples.

`profile` : input, real\*8, dimension (`nprofile`, 2)
:   List of latitude (:,1) and longitude (:,2) coordinates in degrees specifying a single closed curve.

`nprofile` : input, integer
:   The number of coordinates in the curve `profile`.

`np` : input, integer
:   The value of the returned mask at the North pole (90N, 0E). If the North pole is outside of the concentration region, set this to 0; if it is inside the concentration region, set this to 1.

`centralmeridian` : input, optional, integer, default = 0
:   If 1, the curve is assumed to pass through the central meridian: passing from < 360 degrees to > 0 degrees. The curve makes a complete circle about the planet in longitude.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`Curve2Mask` will take a list of latitude and longitude coordinates that represent a single closed curve, and output a mask `mask_dh` that contains 1s and 0s where the grid nodes are inside and outside of the curve, respectively. `mask_dh` must be sampled according to the Driscoll and Healy (1994) sampling theorem with `n` samples in latitude, and either possess `n` samples in longitude (`sampling=1`) or `2n` samples in longitude (`sampling=2`). It is necessary to specify a single point as being inside or outside of the curve, and for this the value at the North pole (90N, 0E) must be specified as either 0 or 1.

This routine saves the three-term recursion factors and square roots of the integers the first time being called. If subsequent calls possess the same value of `lmax`, these will not be recomputed. If you wish to deallocate this memory, which is an array of length `(lmax+1)*(lmax+2)`, recall this routine with `lmax=-1`.

## References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

## See also

[shreturntapersmap](shreturntapersmap.html), [computedmap](computedmap.html)
